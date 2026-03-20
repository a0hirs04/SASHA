#include <cassert>
#include <iostream>
#include <vector>

#include "../Stroma_world/PhysiCell/core/PhysiCell.h"
#include "../Stroma_world/PhysiCell/modules/PhysiCell_standard_modules.h"
#include "../custom_modules/custom.h"

using namespace BioFVM;
using namespace PhysiCell;

int main()
{
    const bool xml_status = load_PhysiCell_config_file("config/PhysiCell_settings.xml");
    assert(xml_status);

    setup_microenvironment();
    create_cell_container_for_microenvironment(microenvironment, 30.0);
    create_cell_types();

    if (parameters.doubles.find_index("tgfb_activation_threshold") >= 0)
    {
        parameters.doubles("tgfb_activation_threshold") = 0.6;
    }
    if (parameters.doubles.find_index("shh_activation_threshold") >= 0)
    {
        parameters.doubles("shh_activation_threshold") = 0.05;
    }

    Cell_Definition* pStroma = find_cell_definition("stromal_cell");
    assert(pStroma != NULL);
    assert(tgfb_index >= 0);
    assert(shh_index >= 0);

    // Rule 22-A: TGF-beta alone above threshold activates.
    Cell* cell_tgfb_only = create_cell(*pStroma);
    cell_tgfb_only->assign_position(std::vector<double>{80.0, 100.0, 0.0});
    const int tgfb_only_acta2 = cell_tgfb_only->custom_data.find_variable_index("acta2_active");
    assert(tgfb_only_acta2 >= 0);
    std::vector<double>& rho_tgfb_only = cell_tgfb_only->nearest_density_vector();
    rho_tgfb_only[tgfb_index] = 0.7;
    rho_tgfb_only[shh_index] = 0.0;
    module3_stromal_activation(cell_tgfb_only, cell_tgfb_only->phenotype, 1.0, ModulePhase::SENSING);
    assert(cell_tgfb_only->custom_data[tgfb_only_acta2] == 1.0);
    std::cout << "PASS Rule22_A_TGFB_alone" << std::endl;

    // Rule 22-B: SHH alone above threshold activates.
    Cell* cell_shh_only = create_cell(*pStroma);
    cell_shh_only->assign_position(std::vector<double>{90.0, 100.0, 0.0});
    const int shh_only_acta2 = cell_shh_only->custom_data.find_variable_index("acta2_active");
    const int shh_only_gli1 = cell_shh_only->custom_data.find_variable_index("gli1_active");
    assert(shh_only_acta2 >= 0);
    assert(shh_only_gli1 >= 0);
    std::vector<double>& rho_shh_only = cell_shh_only->nearest_density_vector();
    rho_shh_only[tgfb_index] = 0.0;
    rho_shh_only[shh_index] = 0.7;
    module3_stromal_activation(cell_shh_only, cell_shh_only->phenotype, 1.0, ModulePhase::SENSING);
    assert(cell_shh_only->custom_data[shh_only_acta2] == 1.0);
    assert(cell_shh_only->custom_data[shh_only_gli1] > 0.9);
    const int shh_only_tgfb = cell_shh_only->custom_data.find_variable_index("tgfb_secretion_active");
    const int shh_only_ecm = cell_shh_only->custom_data.find_variable_index("ecm_production_rate");
    assert(shh_only_tgfb >= 0);
    assert(shh_only_ecm >= 0);
    assert(cell_shh_only->custom_data[shh_only_tgfb] > 0.9);
    assert(cell_shh_only->custom_data[shh_only_ecm] > 0.01);
    std::cout << "PASS Rule22_B_SHH_alone" << std::endl;

    // Rule 22-C: both together below individual level but above combined threshold.
    Cell* cell_combined = create_cell(*pStroma);
    cell_combined->assign_position(std::vector<double>{95.0, 100.0, 0.0});
    const int combined_acta2 = cell_combined->custom_data.find_variable_index("acta2_active");
    assert(combined_acta2 >= 0);
    std::vector<double>& rho_combined = cell_combined->nearest_density_vector();
    rho_combined[tgfb_index] = 0.35;
    rho_combined[shh_index] = 0.30; // each < 0.6, sum = 0.65 > 0.6
    module3_stromal_activation(cell_combined, cell_combined->phenotype, 1.0, ModulePhase::SENSING);
    assert(cell_combined->custom_data[combined_acta2] == 1.0);
    std::cout << "PASS Rule22_C_combined_above_threshold" << std::endl;

    // Rule 22-D: both below combined threshold stay PSC.
    Cell* cell_below = create_cell(*pStroma);
    cell_below->assign_position(std::vector<double>{97.0, 100.0, 0.0});
    const int below_acta2 = cell_below->custom_data.find_variable_index("acta2_active");
    assert(below_acta2 >= 0);
    std::vector<double>& rho_below = cell_below->nearest_density_vector();
    rho_below[tgfb_index] = 0.2;
    rho_below[shh_index] = 0.1; // sum = 0.3 < 0.6
    module3_stromal_activation(cell_below, cell_below->phenotype, 1.0, ModulePhase::SENSING);
    assert(cell_below->custom_data[below_acta2] == 0.0);
    std::cout << "PASS Rule22_D_combined_below_threshold" << std::endl;

    // Case 1: activation + SHH contribution => ACTA2 on (irreversible), GLI1 on.
    Cell* cell_a = create_cell(*pStroma);
    cell_a->assign_position(std::vector<double>{100.0, 100.0, 0.0});

    const int acta2_idx_a = cell_a->custom_data.find_variable_index("acta2_active");
    const int gli1_idx_a = cell_a->custom_data.find_variable_index("gli1_active");
    const int tgfb_sec_idx_a = cell_a->custom_data.find_variable_index("tgfb_secretion_active");
    const int ecm_prod_idx_a = cell_a->custom_data.find_variable_index("ecm_production_rate");
    assert(acta2_idx_a >= 0);
    assert(gli1_idx_a >= 0);
    assert(tgfb_sec_idx_a >= 0);
    assert(ecm_prod_idx_a >= 0);

    std::vector<double>& densities_a = cell_a->nearest_density_vector();
    densities_a[tgfb_index] = 0.5;
    densities_a[shh_index] = 0.3;
    module3_stromal_activation(cell_a, cell_a->phenotype, 1.0, ModulePhase::SENSING);

    assert(cell_a->custom_data[acta2_idx_a] == 1.0);
    assert(cell_a->custom_data[gli1_idx_a] > 0.8);
    assert(cell_a->custom_data[tgfb_sec_idx_a] > 0.0);
    const double ecm_prod_with_shh = cell_a->custom_data[ecm_prod_idx_a];
    assert(ecm_prod_with_shh > 0.0);

    // Case 2: activated CAF with SHH gone => ACTA2 stays on, GLI1 reverts off.
    densities_a[shh_index] = 0.0;
    module3_stromal_activation(cell_a, cell_a->phenotype, 1.0, ModulePhase::SENSING);

    assert(cell_a->custom_data[acta2_idx_a] == 1.0);
    assert(cell_a->custom_data[gli1_idx_a] == 0.0);
    assert(cell_a->custom_data[tgfb_sec_idx_a] == 0.0);
    assert(cell_a->custom_data[ecm_prod_idx_a] < ecm_prod_with_shh);

    // Case 2b: during active SHH inhibition, ACTA2+ CAFs keep partial
    // TGF-beta support even when GLI1 falls, while ECM production remains
    // below the SHH-on state.
    parameters.doubles("shh_inhibition_start_time") = 0.0;
    parameters.doubles("shh_inhibition_strength") = 1.0;
    PhysiCell_globals.current_time = 10.0;
    densities_a[tgfb_index] = 0.0;
    densities_a[shh_index] = 0.0;
    module3_stromal_activation(cell_a, cell_a->phenotype, 1.0, ModulePhase::SENSING);
    assert(cell_a->custom_data[acta2_idx_a] == 1.0);
    assert(cell_a->custom_data[gli1_idx_a] == 0.0);
    assert(cell_a->custom_data[tgfb_sec_idx_a] >= 0.59);
    assert(cell_a->custom_data[ecm_prod_idx_a] < ecm_prod_with_shh);
    parameters.doubles("shh_inhibition_start_time") = 1e18;
    parameters.doubles("shh_inhibition_strength") = 0.0;
    PhysiCell_globals.current_time = 0.0;

    // Case 3: below activation threshold => ACTA2 remains off.
    Cell* cell_b = create_cell(*pStroma);
    cell_b->assign_position(std::vector<double>{140.0, 100.0, 0.0});

    const int acta2_idx_b = cell_b->custom_data.find_variable_index("acta2_active");
    assert(acta2_idx_b >= 0);

    std::vector<double>& densities_b = cell_b->nearest_density_vector();
    densities_b[tgfb_index] = 0.2;
    densities_b[shh_index] = 0.1;
    module3_stromal_activation(cell_b, cell_b->phenotype, 1.0, ModulePhase::SENSING);

    assert(cell_b->custom_data[acta2_idx_b] == 0.0);

    // Rule 23: ACTA2 remains ON after full signal removal over long horizon.
    densities_a[tgfb_index] = 0.0;
    densities_a[shh_index] = 0.0;
    for (int step = 0; step < 200; ++step)
    {
        module3_stromal_activation(cell_a, cell_a->phenotype, 1.0, ModulePhase::SENSING);
    }
    assert(cell_a->custom_data[acta2_idx_a] == 1.0);
    assert(cell_a->custom_data[gli1_idx_a] == 0.0);
    assert(cell_a->custom_data[tgfb_sec_idx_a] == 0.0);
    std::cout << "PASS Rule23_irreversible_ACTA2" << std::endl;

    std::cout << "PASS module3_stromal_activation_test" << std::endl;
    return 0;
}
