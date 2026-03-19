#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "../Stroma_world/PhysiCell/core/PhysiCell.h"
#include "../Stroma_world/PhysiCell/modules/PhysiCell_standard_modules.h"
#include "../custom_modules/custom.h"

using namespace BioFVM;
using namespace PhysiCell;

namespace
{

bool nearly_equal(double a, double b, double eps = 1e-10)
{
    return std::fabs(a - b) <= eps;
}

int voxel_index_for_cell(Cell* pCell)
{
    std::vector<double> pos = pCell->position;
    return microenvironment.nearest_voxel_index(pos);
}

void set_ecm_at_voxel(int voxel_index, double density, double ha_fraction)
{
    microenvironment.density_vector(voxel_index)[ecm_index] = density;
    set_ecm_ha_fraction(voxel_index, ha_fraction);
}

double get_ecm_density_at_voxel(int voxel_index)
{
    return microenvironment.density_vector(voxel_index)[ecm_index];
}

} // namespace

int main()
{
    const bool xml_status = load_PhysiCell_config_file("config/PhysiCell_settings.xml");
    assert(xml_status);

    setup_microenvironment();
    create_cell_container_for_microenvironment(microenvironment, 30.0);
    create_cell_types();

    assert(ecm_index >= 0);

    parameters.doubles("mechanical_compaction_strength") = 0.5;
    parameters.doubles("compaction_ecm_increment") = 0.01;
    parameters.doubles("crowding_base_pressure") = 0.0;

    const double dt = 6.0;

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    assert(pTumor != NULL);

    Cell* cell = create_cell(*pTumor);
    cell->assign_position(std::vector<double>{100.0, 100.0, 0.0});

    const int mech_idx = cell->custom_data.find_variable_index("mechanical_pressure");
    assert(mech_idx >= 0);

    const int voxel = voxel_index_for_cell(cell);
    set_ecm_at_voxel(voxel, 0.6, 0.4);
    cell->state.simple_pressure = 2.0;

    // Test A — Solid stress computed correctly.
    module8_mechanical_compaction(cell, cell->phenotype, dt, ModulePhase::WRITE);
    assert(nearly_equal(cell->custom_data[mech_idx], 0.36));
    std::cout << "PASS Test A" << std::endl;

    // Test B — Compaction increases ECM density.
    assert(nearly_equal(get_ecm_density_at_voxel(voxel), 0.605184));
    std::cout << "PASS Test B" << std::endl;

    // Test C — Compaction does not change HA/collagen ratio.
    assert(nearly_equal(get_ecm_ha_fraction(voxel), 0.4));
    std::cout << "PASS Test C" << std::endl;

    // Test D — No compaction when ECM is absent.
    set_ecm_at_voxel(voxel, 0.0, 0.4);
    cell->state.simple_pressure = 2.0;
    module8_mechanical_compaction(cell, cell->phenotype, dt, ModulePhase::WRITE);
    assert(nearly_equal(get_ecm_density_at_voxel(voxel), 0.0));
    std::cout << "PASS Test D" << std::endl;

    // Test E — ECM clamped at 1.0.
    set_ecm_at_voxel(voxel, 0.99, 0.4);
    cell->state.simple_pressure = 100.0;
    module8_mechanical_compaction(cell, cell->phenotype, dt, ModulePhase::WRITE);
    assert(nearly_equal(get_ecm_density_at_voxel(voxel), 1.0));
    std::cout << "PASS Test E" << std::endl;

    std::cout << "PASS module8_mechanical_compaction_test" << std::endl;
    return 0;
}
