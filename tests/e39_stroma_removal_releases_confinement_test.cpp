#include <cassert>
#include <iostream>
#include <vector>

#include "feedback_loop_common.h"

using namespace BioFVM;
using namespace PhysiCell;
using namespace feedback_loop_common;

int main()
{
    initialize_world();

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    assert(pTumor != NULL);
    assert(ecm_index >= 0);

    // No drug in this rule; isolate mechanical confinement.
    parameters.doubles("base_proliferation_rate") = 0.003;
    parameters.doubles("apoptosis_resistance") = 1.0;
    parameters.doubles("drug_kill_coefficient") = 0.0;
    parameters.doubles("drug_uptake_rate") = 0.0;
    parameters.doubles("efflux_drug_reduction") = 0.0;
    parameters.doubles("go_grow_penalty") = 0.0;
    parameters.doubles("hypoxia_proliferation_modifier") = 0.0;
    parameters.doubles("hypoxia_death_resistance_bonus") = 0.0;

    parameters.doubles("contact_inhibition_threshold") = 0.4;
    parameters.doubles("mechanical_compaction_strength") = 1.0;
    parameters.doubles("compaction_ecm_increment") = 0.0;

    // Keep environment simple and oxygen-rich.
    reset_all_fields(38.0, 0.0, 0.0, 0.0, 0.0, 0.5);

    // Tumor initial condition.
    place_tumor_cluster(pTumor, 20, 0.0, 0.0, 100.0);
    const int tumor_type = pTumor->type;

    // Dense collagen-dominant ring confining expansion.
    set_annulus_ecm(0.0, 0.0, 130.0, 260.0, 0.9, 0.1);

    double t = 0.0;
    const double dt = 6.0;

    const int n0 = static_cast<int>(live_cells_of_type(tumor_type).size());
    advance_steps(40, dt, t);
    const int n100 = static_cast<int>(live_cells_of_type(tumor_type).size());
    const int growth_confined = n100 - n0;

    // Remove ECM barrier completely (stroma removal surrogate).
    const int n_vox = static_cast<int>(microenvironment.number_of_voxels());
    for (int n = 0; n < n_vox; ++n)
    {
        microenvironment.density_vector(n)[ecm_index] = 0.0;
        set_ecm_ha_fraction(n, 0.5);
    }

    advance_steps(40, dt, t);
    const int n200 = static_cast<int>(live_cells_of_type(tumor_type).size());
    const int growth_released = n200 - n100;

    std::cout << "Rule39 metrics"
              << " n0=" << n0
              << " n100_confined=" << n100
              << " n200_after_removal=" << n200
              << " growth_confined=" << growth_confined
              << " growth_after_removal=" << growth_released
              << std::endl;

    assert(growth_released > growth_confined);
    assert(n200 > n100);

    std::cout << "PASS Rule39_stroma_removal_worsens_tumor_growth" << std::endl;
    std::cout << "PASS e39_stroma_removal_releases_confinement_test" << std::endl;
    return 0;
}
