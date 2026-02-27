#include "tumor_cell.h"

#include <algorithm>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

using namespace PhysiCell;

// Substrate indices are resolved once in custom.cpp::setup_microenvironment().
// Keeping them global avoids mismatches between modules.
extern int oxygen_index;
extern int tgfb_index;
extern int shh_index;
extern int ecm_index;
extern int drug_index;

namespace
{

inline double clamp_unit(double x)
{
    if (x < 0.0) return 0.0;
    if (x > 1.0) return 1.0;
    return x;
}

inline double clamp_nonnegative(double x)
{
    return (x < 0.0) ? 0.0 : x;
}

inline double read_density_value(const std::vector<double>& densities, int index)
{
    if (index < 0) return 0.0;
    if (index >= static_cast<int>(densities.size())) return 0.0;
    return densities[index];
}

inline void write_density_value(std::vector<double>& densities, int index, double value)
{
    if (index < 0) return;
    if (index >= static_cast<int>(densities.size())) return;
    densities[index] = value;
}

inline void set_custom_data_if_present(Cell* pCell, const std::string& name, double value)
{
    const int idx = pCell->custom_data.find_variable_index(name);
    if (idx >= 0)
    {
        pCell->custom_data[idx] = value;
    }
}

double get_tumor_base_proliferation_rate(Phenotype& phenotype)
{
    static bool initialized = false;
    static double base_rate = 0.0;

    if (!initialized)
    {
        Cell_Definition* pTumorDef = find_cell_definition("tumor_cell");
        if (pTumorDef != NULL)
        {
            base_rate = pTumorDef->phenotype.cycle.data.transition_rate(0, 0);
        }
        else
        {
            // Fallback: use the first observed phenotype value if the definition lookup fails.
            base_rate = phenotype.cycle.data.transition_rate(0, 0);
        }
        initialized = true;
    }

    return base_rate;
}

double get_tumor_base_apoptosis_rate(Phenotype& phenotype, int apoptosis_index)
{
    static bool initialized = false;
    static double base_rate = 0.0;

    if (!initialized)
    {
        Cell_Definition* pTumorDef = find_cell_definition("tumor_cell");
        if (pTumorDef != NULL)
        {
            const int def_apoptosis_index =
                pTumorDef->phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);
            if (def_apoptosis_index >= 0 &&
                def_apoptosis_index < static_cast<int>(pTumorDef->phenotype.death.rates.size()))
            {
                base_rate = pTumorDef->phenotype.death.rates[def_apoptosis_index];
            }
        }

        if (base_rate <= 0.0 &&
            apoptosis_index >= 0 &&
            apoptosis_index < static_cast<int>(phenotype.death.rates.size()))
        {
            // Secondary fallback if tumor definition is unavailable.
            base_rate = phenotype.death.rates[apoptosis_index];
        }

        initialized = true;
    }

    return base_rate;
}

// Per-cell Boolean network cache keyed by PhysiCell cell ID.
std::unordered_map<int, std::unique_ptr<BooleanNetwork> > g_boolean_networks_by_cell_id;
std::mutex g_boolean_network_mutex;

// Global intervention list, typically populated once at simulation startup.
std::vector<Intervention> g_current_interventions;

// Runtime thresholds loaded from XML user_parameters.
ThresholdConfig g_threshold_config;
bool g_threshold_config_initialized = false;
std::mutex g_threshold_mutex;

} // namespace

BooleanNetwork* get_boolean_network(Cell* pCell, CellType cell_type)
{
    if (pCell == NULL)
    {
        return NULL;
    }

    const int cell_id = pCell->ID;

    std::lock_guard<std::mutex> guard(g_boolean_network_mutex);

    std::unique_ptr<BooleanNetwork>& bn_ptr = g_boolean_networks_by_cell_id[cell_id];
    if (!bn_ptr)
    {
        bn_ptr.reset(new BooleanNetwork());

        // Initialize once per cell with default parameters and cell-type-specific mutant locks.
        GeneParams params;
        bn_ptr->initialize(cell_type, params);

        // Seed network state from current custom_data values for this specific cell.
        bn_ptr->sync_from_cell(pCell);
    }

    return bn_ptr.get();
}

const std::vector<Intervention>& get_current_interventions()
{
    return g_current_interventions;
}

void set_current_interventions(const std::vector<Intervention>& interventions)
{
    // Intended to be called during setup before parallel simulation starts.
    g_current_interventions = interventions;
}

const ThresholdConfig& get_threshold_config()
{
    std::lock_guard<std::mutex> guard(g_threshold_mutex);

    if (!g_threshold_config_initialized)
    {
        g_threshold_config = ThresholdConfig::load_from_xml();
        g_threshold_config_initialized = true;
    }

    return g_threshold_config;
}

void set_threshold_config(const ThresholdConfig& cfg)
{
    std::lock_guard<std::mutex> guard(g_threshold_mutex);
    g_threshold_config = cfg;
    g_threshold_config_initialized = true;
}

void tumor_phenotype_update(Cell* pCell, Phenotype& phenotype, double dt)
{
    // Defensive guards: null cells, dead cells, or non-positive dt should not be updated.
    if (pCell == NULL) return;
    if (phenotype.death.dead) return;
    if (dt <= 0.0) return;

    // -------------------------------------------------------------------------
    // STEP 1. Read microenvironment values at the cell's current voxel.
    // -------------------------------------------------------------------------
    std::vector<double>& densities = pCell->nearest_density_vector();

    const double oxygen        = read_density_value(densities, oxygen_index);
    const double tgfb          = read_density_value(densities, tgfb_index);
    const double shh           = read_density_value(densities, shh_index);
    const double drug          = read_density_value(densities, drug_index);
    // ECM density read once here; used for physical barrier effects on
    // proliferation, motility, and drug delivery throughout this function.
    const double local_ecm_val = read_density_value(densities, ecm_index);

    // -------------------------------------------------------------------------
    // STEP 2. Update the per-cell Boolean network state.
    // -------------------------------------------------------------------------
    BooleanNetwork* bn = get_boolean_network(pCell, CellType::TUMOR);
    if (bn == NULL) return;

    // Keep the network synchronized with cell-level custom_data state.
    bn->sync_from_cell(pCell);

    // Use runtime XML-backed thresholds (single source of truth) during update.
    const ThresholdConfig& cfg = get_threshold_config();
    bn->update(dt, oxygen, tgfb, shh, drug, cfg);

    // Apply globally loaded interventions after intrinsic network dynamics.
    bn->apply_interventions(get_current_interventions());

    // Persist updated gene states back to PhysiCell custom_data.
    bn->sync_to_cell(pCell);

    // -------------------------------------------------------------------------
    // STEP 3. Dominance voting: compute qualitative axis outcomes.
    //         (Gene State Representation & Combination Rules v1.0, Part 3-4)
    // -------------------------------------------------------------------------

    // Short aliases (continuous gene activities after ODE update).
    const double zeb1    = clamp_unit((*bn)[ZEB1]);
    const double cdh1    = clamp_unit((*bn)[CDH1]);
    const double snai1   = clamp_unit((*bn)[SNAI1]);
    const double mmp2    = clamp_unit((*bn)[MMP2]);
    const double tgfb1   = clamp_unit((*bn)[TGFB1]);
    const double shh_expr= clamp_unit((*bn)[SHH]);
    const double nrf2    = clamp_unit((*bn)[NRF2]);
    const double abcb1   = clamp_unit((*bn)[ABCB1]);
    const double bcl_xl  = clamp_unit((*bn)[BCL_XL]);
    const double cdkn2a  = clamp_unit((*bn)[CDKN2A]);
    const double tp53    = clamp_unit((*bn)[TP53]);
    const double myc     = clamp_unit((*bn)[MYC]);

    // Compute qualitative axis outcomes via dominance voting.
    // compute_axis_outcome() implements the signal contribution table (Part 3)
    // and the dominance voting rules (Part 4, Steps 1–2) from the document.
    AxisOutcome growth_outcome   = bn->compute_axis_outcome(FunctionalAxis::GROWTH,
                                                            tgfb, drug);
    AxisOutcome death_outcome    = bn->compute_axis_outcome(FunctionalAxis::DEATH,
                                                            tgfb, drug);
    AxisOutcome braking_outcome  = bn->compute_axis_outcome(FunctionalAxis::BRAKING,
                                                            tgfb, drug);
    AxisOutcome invasion_outcome = bn->compute_axis_outcome(FunctionalAxis::INVASION,
                                                            tgfb, drug);

    // -------------------------------------------------------------------------
    // STEP 4. Apply veto rules (Part 4, Step 3 — hard logical overrides).
    //         Veto rules enforce biological constraints that override voting.
    // -------------------------------------------------------------------------

    // Veto 1: CDKN2A=- AND TP53=- AND MYC=+
    //   → Cell Cycle Braking = NONE (arrest is impossible regardless of signals).
    //   All three major brakes destroyed + override active.
    //   This is the canonical PDAC state. PhysiCell: proliferation is unconstrained
    //   by intracellular signals (only mechanical limits apply).
    if (cdkn2a < 0.1 && tp53 < 0.1 && myc > 0.75)
    {
        braking_outcome = AxisOutcome::NONE;
    }

    // Veto 2: BCL-XL=+ AND BAX=WT (functionally outcompeted by BCL-XL)
    //   → Death Resistance floor = HIGH.
    //   Anti-apoptotic dominant. Apoptosis requires external intervention.
    //   BAX is not in the 21-gene panel; BCL-XL activity serves as the gate.
    if (bcl_xl > 0.75 && static_cast<int>(death_outcome) < static_cast<int>(AxisOutcome::HIGH))
    {
        death_outcome = AxisOutcome::HIGH;
    }

    // Veto 3: SMAD4=- AND TGFB1 present → TGF-β brake on tumor = DISABLED.
    //   Already structurally enforced in compute_tumor_targets() (RB1 rule):
    //   SMAD4=0 makes the SMAD4*tgfb product identically zero → RB1 cannot rise.
    //   No additional override needed here; braking_outcome already reflects this.

    // Veto 4: ZEB1=active AND CDH1=repressed → Motility floor = HIGH.
    //   Full EMT. Cell is mesenchymal regardless of other adhesion signals.
    //   Clinical note: this is the invasive front phenotype — the most drug-resistant zone.
    if (zeb1 > 0.5 && cdh1 < 0.3 &&
        static_cast<int>(invasion_outcome) < static_cast<int>(AxisOutcome::HIGH))
    {
        invasion_outcome = AxisOutcome::HIGH;
    }

    // Veto 5: ABCB1=active AND drug absent → ABCB1 resistance = ZERO.
    //   Efflux pump has nothing to pump. (Handled in drug_sensitivity_stored below.)

    // -------------------------------------------------------------------------
    // STEP 5. Map axis outcomes to PhysiCell phenotype parameters via bins.
    //         Numbers enter ONLY here. Bins are literature-backed (PhenoParamBins).
    //         Physical modifiers (ECM, solid stress) are applied as multipliers.
    // -------------------------------------------------------------------------

    // ---- 5a) PROLIFERATION --------------------------------------------------
    // Growth bin: converts GROWTH axis outcome to a base proliferation rate.
    const double growth_bin = PhenoParamBins::get(PhenoParamBins::PROLIF, growth_outcome);

    // Braking multiplier: converts BRAKING outcome to a growth suppressor.
    // NONE braking (veto 1) → multiplier = 1.0 (unstoppable).
    // VERY_HIGH braking → multiplier ≈ 0 (fully arrested).
    // Linear: braking_mult = 1 - (braking_int / 6).
    const int braking_int = static_cast<int>(braking_outcome);
    const double braking_mult = (braking_outcome == AxisOutcome::NONE)
        ? 1.0
        : clamp_unit(1.0 - braking_int * (1.0 / 6.0));

    // Contact inhibition: physical crowding suppresses proliferation in brake-
    // intact cells. When braking_outcome=NONE (veto 1 fired), CDKN2A-/TP53-
    // cells bypass crowding signals via constitutive growth factor drive.
    const double neighbor_fraction =
        clamp_unit(static_cast<double>(pCell->state.neighbors.size()) / 8.0);
    double contact_mult = 1.0;
    if (braking_outcome != AxisOutcome::NONE &&
        neighbor_fraction > cfg.contact_inhibition_threshold)
    {
        contact_mult = 0.1;
    }

    // ECM physical constraint (vismodegib failure mechanism — Phase 1.2 Caveat #6).
    // Dense desmoplastic stroma compresses tumor mass, raising interstitial fluid
    // pressure and limiting outward expansion. Removing this containment without
    // a cytotoxic agent → tumor expansion into freed space → worse outcome.
    // Clinical reference: Rhim et al. 2014, doi:10.1016/j.ccr.2014.04.009.
    // At ECM=0.1 (baseline): ~2.5% constraint. At ECM=0.8 (desmoplastic): ~20%.
    const double ecm_brake = clamp_unit(1.0 - 0.25 * local_ecm_val);

    // ---- TRAIT 8: Passive compaction / solid stress -------------------------
    // Mechanical compression from tumor cell growth against stiff desmoplastic
    // ECM. This is DISTINCT from the biochemical ECM delivery barrier (5c) and
    // the ECM physical constraint above — solid stress is the MECHANICAL load
    // from cells that cannot expand due to matrix resistance.
    //
    // pCell->state.simple_pressure: normalized crowding pressure from neighbor
    // overlap. Dense ECM (high collagen/HA stiffness) resists strain relaxation,
    // amplifying and retaining solid stress rather than allowing dissipation.
    //
    // In PDAC, solid stress from growth against desmoplastic stroma accounts for
    // ~10-20% of total growth limitation and compresses tumor vasculature
    // independently of interstitial fluid pressure.
    // Reference: Stylianopoulos et al. 2012, PNAS doi:10.1073/pnas.1213353109
    // Without this term, the barrier is softer than in vivo and the EA will find
    // drug-delivery solutions that would not work clinically.
    const double simple_pressure = clamp_nonnegative(pCell->state.simple_pressure);
    // ECM stiffness amplifier: dense collagen/HA resists deformation, retaining
    // mechanical stress rather than allowing viscoelastic dissipation.
    // At ECM=0.0: solid_stress = simple_pressure (cells only).
    // At ECM=0.8: solid_stress = 1.8 × simple_pressure (stiff matrix contribution).
    const double solid_stress = simple_pressure * (1.0 + local_ecm_val);
    // Saturating brake: 35% max suppression at high solid stress (cells adapt
    // via cytoskeletal remodeling, preventing complete growth arrest from
    // mechanical stress alone — unlike chemical veto rule 1).
    const double solid_stress_brake = clamp_unit(1.0 - 0.35 * std::min(1.0, solid_stress));

    // Final proliferation rate: growth_bin × braking × contact × ECM × solid_stress.
    const double final_prolif = growth_bin * braking_mult * contact_mult
                              * ecm_brake * solid_stress_brake;
    phenotype.cycle.data.transition_rate(0, 0) = std::max(0.0, final_prolif);

    // ---- 5b) APOPTOSIS RATE -------------------------------------------------
    // Death outcome represents RESISTANCE; invert to get apoptosis rate.
    // VERY_HIGH resistance → NONE apoptosis; NONE resistance → VERY_HIGH apoptosis.
    const int apoptosis_index =
        phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);
    if (apoptosis_index >= 0 &&
        apoptosis_index < static_cast<int>(phenotype.death.rates.size()))
    {
        // Invert: resistance 6→apoptosis 0, resistance 0→apoptosis 6.
        int death_int = static_cast<int>(death_outcome);
        int apopt_int = 6 - death_int;
        if (apopt_int < 0) apopt_int = 0;
        if (apopt_int > 6) apopt_int = 6;
        const double base_apoptosis_rate =
            PhenoParamBins::get(PhenoParamBins::APOPTOSIS,
                                static_cast<AxisOutcome>(apopt_int));
        phenotype.death.rates[apoptosis_index] = base_apoptosis_rate;

        // Drug-induced death (additive, on top of basal resistance).
        // Veto 5: ABCB1 counts toward resistance ONLY when drug is present.
        if (drug > 0.01)
        {
            // Drug sensitivity (NRF2 + ABCB1 resistance, both induced by drug).
            // Clamped to [0.05, 1.0]: floor prevents complete drug immunity.
            double drug_sensitivity = 1.0 - 0.4 * nrf2 - 0.4 * abcb1;
            drug_sensitivity = std::max(0.05, std::min(1.0, drug_sensitivity));

            // ECM delivery barrier (PEGPH20 failure mechanism — Phase 1.2 Caveat):
            // Dense HA + collagen raises IFP, reducing convective drug transport.
            //   (a) PEGPH20 (HA depletion) improves gemcitabine in preclinical PDAC.
            //   (b) PEGPH20 alone fails clinically — COL1A1 barrier remains intact.
            // At ECM=0.1: 5% barrier. At ECM=0.8: 40% barrier.
            const double ecm_delivery_barrier = clamp_unit(1.0 - 0.5 * local_ecm_val);
            const double effective_drug = clamp_nonnegative(drug) * drug_sensitivity
                                        * ecm_delivery_barrier;
            phenotype.death.rates[apoptosis_index] += effective_drug * 0.001;
        }
    }

    // ---- 5c) MOTILITY / EMT -------------------------------------------------
    // TRAIT 5: Inducible spatial EMT — only hypoxic cells go mesenchymal.
    // This is spatial by design: oxygen varies across the domain (PhysiCell PDE).
    // HIF1A activates ZEB1 only when oxygen < HYPOXIA_HIF1A_THRESHOLD.
    // Result: normoxic periphery stays epithelial; hypoxic core goes mesenchymal.
    // The invasive front is the most immune-resistant zone — exactly the
    // layered defense real PDAC presents. No additional code needed here;
    // the spatial oxygen field drives the HIF1A → ZEB1 → CDH1 chain automatically.
    //
    // SNAI1 can independently trigger mesenchymal motility before ZEB1 rises
    // (fast responder to TGF-β: hours vs. days for ZEB1).
    const bool is_mesenchymal = (zeb1  > cfg.emt_zeb1_threshold ||
                                 snai1 > cfg.emt_zeb1_threshold);
    set_custom_data_if_present(pCell, "is_mesenchymal", is_mesenchymal ? 1.0 : 0.0);

    // Migration speed from invasion axis bin + ECM steric factor.
    // ECM steric barrier: dense HA and collagen impede cell movement via
    // increased matrix viscosity and entanglement. Mesenchymal cells use MMP2
    // to carve invasion paths, but this is rate-limited → modeled as a speed cap.
    const double speed_bin = PhenoParamBins::get(PhenoParamBins::MIGRATION,
                                                  invasion_outcome);
    const double ecm_motility_factor = clamp_unit(1.0 - 0.5 * local_ecm_val);
    // floor at 0.01 µm/min: cells always retain minimal basal motility
    phenotype.motility.migration_speed =
        std::max(0.01, speed_bin * ecm_motility_factor);

    // CDH1 (E-cadherin) → cell-cell adhesion: graded mapping captures partial EMT.
    // CDH1=1.0 (epithelial): adhesion = 0.4 (strong, contact-inhibited).
    // CDH1=0.0 (mesenchymal): adhesion = 0.1 (weak, migratory).
    // Partial EMT states (CDH1 ≈ 0.5) are common in PDAC; intermediate values.
    phenotype.mechanics.cell_cell_adhesion_strength = 0.1 + 0.3 * cdh1;

    // ---- 5d) ECM DEGRADATION (local voxel update) ---------------------------
    if (ecm_index >= 0)
    {
        const double local_ecm = read_density_value(densities, ecm_index);
        const double degraded_ecm = clamp_unit(local_ecm - mmp2 * 0.01 * dt);
        write_density_value(densities, ecm_index, degraded_ecm);
    }

    // ---- 5e) SECRETION -------------------------------------------------------
    // Expression proxies (gene activities) map to secretion rates.
    // TRAIT 1: KRAS constitutively drives TGF-β and SHH secretion (KRAS=1, locked).
    // This is the primary mechanism for stroma self-assembly (Example 4 in the doc).
    const double tgfb_secretion = clamp_nonnegative(tgfb1 * 0.05);
    const double shh_secretion  = clamp_nonnegative(shh_expr * 0.02);

    if (tgfb_index >= 0 &&
        tgfb_index < static_cast<int>(phenotype.secretion.secretion_rates.size()))
    {
        phenotype.secretion.secretion_rates[tgfb_index] = tgfb_secretion;
    }
    if (shh_index >= 0 &&
        shh_index < static_cast<int>(phenotype.secretion.secretion_rates.size()))
    {
        phenotype.secretion.secretion_rates[shh_index] = shh_secretion;
    }

    // ---- 5f) DRUG SENSITIVITY (Veto 5 implementation) -----------------------
    // ABCB1 is an INDUCED resistance gene: efflux pump has nothing to pump when
    // no cytotoxic substrate is present. Storing 1.0 (fully sensitive) when drug
    // is absent prevents the EA from treating ABCB1 as a constitutive survival
    // gene and wasting generations targeting it before treatment begins.
    // (Gene State Representation, Veto Rule 5: ABCB1=active AND drug absent → ZERO)
    const double drug_sensitivity_stored = (drug > 0.01)
        ? std::max(0.05, std::min(1.0, 1.0 - 0.4 * nrf2 - 0.4 * abcb1))
        : 1.0;
    set_custom_data_if_present(pCell, "drug_sensitivity", drug_sensitivity_stored);

    // Store axis outcomes for EA fitness evaluation and behavioral test validation.
    set_custom_data_if_present(pCell, "axis_growth",   static_cast<double>(growth_outcome));
    set_custom_data_if_present(pCell, "axis_death",    static_cast<double>(death_outcome));
    set_custom_data_if_present(pCell, "axis_braking",  static_cast<double>(braking_outcome));
    set_custom_data_if_present(pCell, "axis_invasion", static_cast<double>(invasion_outcome));
    set_custom_data_if_present(pCell, "solid_stress",  solid_stress);
}
