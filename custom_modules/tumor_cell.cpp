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
    // STEP 3. Map gene states to PhysiCell phenotype parameters.
    // -------------------------------------------------------------------------

    // Short aliases for readability; clamp each to [0,1] to avoid propagation
    // of any numerical drift into phenotype rules.
    const double kras    = clamp_unit((*bn)[KRAS]);
    const double myc     = clamp_unit((*bn)[MYC]);
    const double ccnd1   = clamp_unit((*bn)[CCND1]);
    const double rb1     = clamp_unit((*bn)[RB1]);
    const double cdkn2a  = clamp_unit((*bn)[CDKN2A]);
    const double smad4   = clamp_unit((*bn)[SMAD4]);
    const double bcl_xl  = clamp_unit((*bn)[BCL_XL]);
    const double snai1   = clamp_unit((*bn)[SNAI1]);
    const double tp53    = clamp_unit((*bn)[TP53]);
    const double zeb1    = clamp_unit((*bn)[ZEB1]);
    const double mmp2    = clamp_unit((*bn)[MMP2]);
    const double tgfb1   = clamp_unit((*bn)[TGFB1]);
    const double shh_expr= clamp_unit((*bn)[SHH]);
    const double nrf2    = clamp_unit((*bn)[NRF2]);
    const double abcb1   = clamp_unit((*bn)[ABCB1]);

    // 3a) PROLIFERATION RATE
    // Growth drive: KRAS sets the base signal, MYC amplifies it, CCND1 gates
    // the cell cycle at G1/S. In canonical PDAC (KRAS=1, MYC≈1, CCND1≈1,
    // CDKN2A=0), growth_drive approaches 1.0.
    const double growth_drive = clamp_unit(kras * (0.4 + 0.35 * myc + 0.25 * ccnd1));

    // RB1 is a brake on cell-cycle progression. High RB1 lowers effective growth.
    double brake = clamp_unit(1.0 - 0.8 * rb1);

    // Approximate contact inhibition from local neighbor count.
    const double neighbor_fraction =
        clamp_unit(static_cast<double>(pCell->state.neighbors.size()) / 8.0);

    if (cdkn2a < 0.1 && smad4 < 0.1)
    {
        // Both brakes genetically lost: enforce only a weak minimum brake.
        brake = std::max(brake, 0.2);
    }
    else
    {
        // Functional brake circuitry: crowding strongly suppresses proliferation.
        if (neighbor_fraction > cfg.contact_inhibition_threshold)
        {
            brake *= 0.1;
        }
    }

    // Go-vs-Grow tradeoff: ZEB1-high (mesenchymal) cells divert transcriptional
    // resources from the proliferative program to the motility/invasion program.
    // ZEB1 directly represses several MYC-target cell-cycle genes (CCND1, CDK4)
    // and disrupts mitotic spindle assembly in the high-ZEB1 state.
    // Ceiling at 35% reduction: partial EMT cells maintain moderate proliferation
    // alongside motility; full EMT (ZEB1=1) trades ~35% growth for peak invasion.
    if (zeb1 > 0.3)
    {
        brake *= clamp_unit(1.0 - 0.35 * zeb1);
    }

    // ECM physical constraint: dense desmoplastic stroma mechanically compresses
    // the tumor mass, elevating interstitial fluid pressure and limiting outward
    // expansion. Releasing this constraint (e.g., SHH inhibition → CAF depletion
    // → ECM regression) paradoxically accelerates tumor growth by removing the
    // physical containment — the vismodegib failure mechanism (Caveat #6).
    // Clinical reference: Rhim et al. 2014, doi:10.1016/j.ccr.2014.04.009.
    // At baseline ECM=0.1: ~2.5% brake. At desmoplastic ECM=0.8: ~20% brake.
    brake *= clamp_unit(1.0 - 0.25 * local_ecm_val);

    const double base_prolif_rate = get_tumor_base_proliferation_rate(phenotype);
    const double mapped_prolif_rate =
        clamp_nonnegative(base_prolif_rate * growth_drive * brake);
    phenotype.cycle.data.transition_rate(0, 0) = mapped_prolif_rate;

    // 3b) APOPTOSIS RATE
    const int apoptosis_index =
        phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);
    if (apoptosis_index >= 0 &&
        apoptosis_index < static_cast<int>(phenotype.death.rates.size()))
    {
        // BCL-XL suppresses apoptosis; TP53 loss (=0) blunts the death signal.
        // death_suppression: when BCL_XL=0.8 (PDAC default) → 0.2; when
        // BCL_XL is therapeutically reduced to 0 → 1.0 (full apoptosis).
        const double death_suppression = clamp_unit(1.0 - bcl_xl);

        // p53_factor: TP53 loss sets a low floor (0.2); therapeutic TP53
        // restoration (tp53→1.0) restores full apoptotic signaling.
        const double p53_factor = clamp_unit(0.2 + 0.8 * tp53);

        const double base_apoptosis_rate =
            get_tumor_base_apoptosis_rate(phenotype, apoptosis_index);

        double mapped_apoptosis_rate =
            clamp_nonnegative(base_apoptosis_rate * death_suppression * p53_factor);
        phenotype.death.rates[apoptosis_index] = mapped_apoptosis_rate;

        // 3g) DRUG-INDUCED DEATH
        // Drug kill is added only when local drug concentration is meaningful.
        if (drug > 0.01)
        {
            // Drug sensitivity clamped to [0.05, 1.0] to avoid complete immunity.
            double drug_sensitivity = 1.0 - 0.4 * nrf2 - 0.4 * abcb1;
            drug_sensitivity = std::max(0.05, std::min(1.0, drug_sensitivity));

            // ECM delivery barrier: dense HA (stromal HAS2) and collagen (COL1A1)
            // compress tumor vasculature, raising interstitial fluid pressure (IFP)
            // and reducing convective drug transport to tumor cells.
            // This explains two clinical observations:
            //   (a) PEGPH20 (HA depletion) improves gemcitabine in preclinical PDAC.
            //   (b) PEGPH20 alone fails clinically — the COL1A1 collagen barrier
            //       remains intact, so drug delivery improves only partially.
            // At ECM=0.1 (baseline): 5% barrier. At ECM=0.8 (desmoplastic): 40% barrier.
            const double ecm_delivery_barrier = clamp_unit(1.0 - 0.5 * local_ecm_val);
            const double effective_drug = clamp_nonnegative(drug) * drug_sensitivity
                                        * ecm_delivery_barrier;
            phenotype.death.rates[apoptosis_index] += effective_drug * 0.001;
        }
    }

    // 3c) MOTILITY (EMT state via ZEB1 or SNAI1)
    // SNAI1 can independently trigger mesenchymal motility before ZEB1 rises.
    const bool is_mesenchymal = (zeb1  > cfg.emt_zeb1_threshold ||
                                 snai1 > cfg.emt_zeb1_threshold);
    set_custom_data_if_present(pCell, "is_mesenchymal", is_mesenchymal ? 1.0 : 0.0);

    // ECM steric barrier reduces migration speed in both EMT states.
    // Dense HA (HAS2) and collagen (COL1A1) physically impede cell movement via
    // increased matrix viscosity and entanglement. Mesenchymal cells use MMP2 to
    // carve invasion paths but this is rate-limited — modeled as a speed cap.
    // At ECM=0.1 (baseline): 5% speed reduction. At ECM=0.8: 40% reduction.
    const double ecm_motility_factor = clamp_unit(1.0 - 0.5 * local_ecm_val);
    if (is_mesenchymal)
    {
        // floor at 0.1 µm/min: mesenchymal cells can still invade even dense ECM
        phenotype.motility.migration_speed =
            std::max(0.1, 1.0 * ecm_motility_factor);
        phenotype.mechanics.cell_cell_adhesion_strength = 0.1;    // weaker adhesion
    }
    else
    {
        // floor at 0.05 µm/min: epithelial cells barely move in thick ECM
        phenotype.motility.migration_speed =
            std::max(0.05, 0.25 * ecm_motility_factor);
        phenotype.mechanics.cell_cell_adhesion_strength = 0.4;    // stronger adhesion
    }

    // 3d) ECM DEGRADATION (local voxel update)
    if (ecm_index >= 0)
    {
        const double local_ecm = read_density_value(densities, ecm_index);
        const double degraded_ecm = clamp_unit(local_ecm - mmp2 * 0.01 * dt);
        write_density_value(densities, ecm_index, degraded_ecm);
    }

    // 3e) SECRETION
    // Expression proxies map to secretion rates for paracrine signaling fields.
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

    // 3f) DRUG SENSITIVITY — stored conditionally on drug presence (Caveat #4).
    // ABCB1 is an INDUCED resistance gene: the efflux pump only matters when
    // cytotoxic drug substrate is present to export. Storing 1.0 (fully sensitive)
    // when no drug is present prevents the EA from treating ABCB1 as a constitutive
    // survival factor and wasting generations targeting it before treatment begins.
    // When drug IS present, the stored value correctly reflects NRF2 + ABCB1
    // resistance for fitness evaluation and downstream analysis.
    const double drug_sensitivity_stored = (drug > 0.01)
        ? std::max(0.05, std::min(1.0, 1.0 - 0.4 * nrf2 - 0.4 * abcb1))
        : 1.0;
    set_custom_data_if_present(pCell, "drug_sensitivity", drug_sensitivity_stored);
}
