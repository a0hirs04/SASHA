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

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    assert(pTumor != NULL);

    Cell* pCell = create_cell(*pTumor);
    pCell->assign_position(std::vector<double>{100.0, 100.0, 0.0});

    assert(oxygen_index >= 0);
    const int hif1a_active_index = pCell->custom_data.find_variable_index("hif1a_active");
    assert(hif1a_active_index >= 0);

    std::vector<double>& densities = pCell->nearest_density_vector();

    densities[oxygen_index] = 0.02;
    module1_oxygen_sensing(pCell, pCell->phenotype, 1.0, ModulePhase::SENSING);
    assert(pCell->custom_data[hif1a_active_index] == 1.0);

    densities[oxygen_index] = 0.08;
    module1_oxygen_sensing(pCell, pCell->phenotype, 1.0, ModulePhase::SENSING);
    assert(pCell->custom_data[hif1a_active_index] == 0.0);

    std::cout << "PASS module1_oxygen_sensing_test" << std::endl;
    return 0;
}
