#include "stromal_cell.h"

#include "tumor_cell.h"

#include <algorithm>
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

double get_stroma_base_proliferation_rate(Phenotype& phenotype)
{
    static bool initialized = false;
    static double base_rate = 0.0;

    if (!initialized)
    {
        Cell_Definition* pStromaDef = find_cell_definition("stromal_cell");
        if (pStromaDef != NULL)
        {
            base_rate = pStromaDef->phenotype.cycle.data.transition_rate(0, 0);
        }
        else
        {
            // Fallback to the first observed cell-level value if definition lookup fails.
            base_rate = phenotype.cycle.data.transition_rate(0, 0);
        }
        initialized = true;
    }

    return base_rate;
}

} // namespace

void stromal_phenotype_update(Cell* pCell, Phenotype& phenotype, double dt)
{
    // Defensive guards.
    if (pCell == NULL) return;
    if (phenotype.death.dead) return;
    if (dt <= 0.0) return;

    // -------------------------------------------------------------------------
    // STEP 1. Read local microenvironment signals.
    // -------------------------------------------------------------------------
    std::vector<double>& densities = pCell->nearest_density_vector();

    const double oxygen = read_density_value(densities, oxygen_index);
    const double tgfb   = read_density_value(densities, tgfb_index);
    const double shh    = read_density_value(densities, shh_index);
    const double drug   = read_density_value(densities, drug_index);

    // -------------------------------------------------------------------------
    // STEP 2. Update per-cell Boolean network state.
    // -------------------------------------------------------------------------
    BooleanNetwork* bn = get_boolean_network(pCell, CellType::STROMA);
    if (bn == NULL) return;

    bn->sync_from_cell(pCell);
    bn->update(dt, oxygen, tgfb, shh, drug, get_threshold_config());
    bn->apply_interventions(get_current_interventions());
    bn->sync_to_cell(pCell);

    // -------------------------------------------------------------------------
    // STEP 3. Map gene states to stromal phenotype behavior.
    // -------------------------------------------------------------------------
    const double acta2  = clamp_unit((*bn)[ACTA2]);
    const double gli1   = clamp_unit((*bn)[GLI1]);
    const double has2   = clamp_unit((*bn)[HAS2]);
    const double col1a1 = clamp_unit((*bn)[COL1A1]);
    const double tgfb1  = clamp_unit((*bn)[TGFB1]);
    const double mmp2   = clamp_unit((*bn)[MMP2]);

    // 3a) ACTIVATION STATE
    const bool is_activated = (acta2 > 0.5);
    set_custom_data_if_present(pCell, "is_activated", is_activated ? 1.0 : 0.0);

    // 3b) PROLIFERATION
    // Quiescent PSCs divide slowly; activated CAFs proliferate faster.
    const double activation_prolif = 0.1 + 0.9 * acta2;
    const double hedgehog_boost = 1.0 + 0.5 * gli1;
    const double base_stroma_prolif_rate = get_stroma_base_proliferation_rate(phenotype);
    phenotype.cycle.data.transition_rate(0, 0) =
        clamp_nonnegative(base_stroma_prolif_rate * activation_prolif * hedgehog_boost);

    // 3c) ECM PRODUCTION
    // ACTA2-gated ECM production from HAS2 (HA) and COL1A1 (collagen).
    const double ecm_production = (0.5 * has2 + 0.5 * col1a1) * acta2;
    set_custom_data_if_present(pCell, "ecm_production_rate", ecm_production);

    // Update local ECM field with production first.
    double local_ecm = read_density_value(densities, ecm_index);
    if (ecm_index >= 0)
    {
        local_ecm = std::min(1.0, local_ecm + ecm_production * 0.005 * dt);
        write_density_value(densities, ecm_index, local_ecm);
    }

    // 3d) MOTILITY
    // Activated CAFs are more contractile and motile than quiescent PSCs.
    if (is_activated)
    {
        phenotype.motility.migration_speed = 0.3;
    }
    else
    {
        phenotype.motility.migration_speed = 0.05;
    }

    // 3e) SECRETION (paracrine/autocrine TGF-beta feedback).
    const double tgfb_secretion = clamp_nonnegative(tgfb1 * 0.03);
    if (tgfb_index >= 0 &&
        tgfb_index < static_cast<int>(phenotype.secretion.secretion_rates.size()))
    {
        phenotype.secretion.secretion_rates[tgfb_index] = tgfb_secretion;
    }

    // 3g) DRUG DIFFUSION MODULATION SUPPORT
    // Store local ECM as an optional custom output for fitness/analysis pipelines.
    set_custom_data_if_present(pCell, "local_ecm_density", local_ecm);
}
