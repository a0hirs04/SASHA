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

    // Case 1: activation + SHH contribution => ACTA2 on (irreversible), GLI1 on.
    Cell* cell_a = create_cell(*pStroma);
    cell_a->assign_position(std::vector<double>{100.0, 100.0, 0.0});

    const int acta2_idx_a = cell_a->custom_data.find_variable_index("acta2_active");
    const int gli1_idx_a = cell_a->custom_data.find_variable_index("gli1_active");
    assert(acta2_idx_a >= 0);
    assert(gli1_idx_a >= 0);

    std::vector<double>& densities_a = cell_a->nearest_density_vector();
    densities_a[tgfb_index] = 0.5;
    densities_a[shh_index] = 0.3;
    module3_stromal_activation(cell_a, cell_a->phenotype, 1.0, ModulePhase::SENSING);

    assert(cell_a->custom_data[acta2_idx_a] == 1.0);
    assert(cell_a->custom_data[gli1_idx_a] == 1.0);

    // Case 2: activated CAF with SHH gone => ACTA2 stays on, GLI1 reverts off.
    densities_a[shh_index] = 0.0;
    module3_stromal_activation(cell_a, cell_a->phenotype, 1.0, ModulePhase::SENSING);

    assert(cell_a->custom_data[acta2_idx_a] == 1.0);
    assert(cell_a->custom_data[gli1_idx_a] == 0.0);

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

    std::cout << "PASS module3_stromal_activation_test" << std::endl;
    return 0;
}
