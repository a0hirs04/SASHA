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

    parameters.doubles("ecm_production_rate_base") = 0.01;
    parameters.doubles("ecm_production_rate_boosted") = 0.015;
    parameters.doubles("ecm_ha_fraction_default") = 0.6;
    parameters.doubles("mmp2_degradation_rate") = 0.005;
    parameters.doubles("ecm_natural_decay_rate") = 0.0001;

    const double dt = 6.0;

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    Cell_Definition* pStroma = find_cell_definition("stromal_cell");
    assert(pTumor != NULL);
    assert(pStroma != NULL);

    intervention_state.ha_degrade_active = false;
    intervention_state.col_degrade_active = false;
    intervention_state.ha_degrade_strength = 0.0;
    intervention_state.col_degrade_strength = 0.0;

    // Test A — Activated CAF produces ECM (no GLI1).
    Cell* caf_a = create_cell(*pStroma);
    caf_a->assign_position(std::vector<double>{100.0, 100.0, 0.0});
    const int a_acta2 = caf_a->custom_data.find_variable_index("acta2_active");
    const int a_gli1 = caf_a->custom_data.find_variable_index("gli1_active");
    assert(a_acta2 >= 0);
    assert(a_gli1 >= 0);
    caf_a->custom_data[a_acta2] = 1.0;
    caf_a->custom_data[a_gli1] = 0.0;
    const int voxel_a = voxel_index_for_cell(caf_a);
    set_ecm_at_voxel(voxel_a, 0.0, 0.5);
    module6_ecm_production(caf_a, caf_a->phenotype, dt, ModulePhase::WRITE);
    assert(nearly_equal(get_ecm_density_at_voxel(voxel_a), 0.06));
    assert(nearly_equal(get_ecm_ha_fraction(voxel_a), 0.6));
    std::cout << "PASS Test A" << std::endl;
    std::cout << "PASS Rule25_caf_produces_HA_and_collagen" << std::endl;

    // Test B — GLI1 boosts production.
    caf_a->custom_data[a_gli1] = 1.0;
    set_ecm_at_voxel(voxel_a, 0.0, 0.5);
    module6_ecm_production(caf_a, caf_a->phenotype, dt, ModulePhase::WRITE);
    assert(nearly_equal(get_ecm_density_at_voxel(voxel_a), 0.09));
    assert(get_ecm_density_at_voxel(voxel_a) > 0.06);
    std::cout << "PASS Test B" << std::endl;
    std::cout << "PASS Rule26_gli1_boosts_ecm_production" << std::endl;

    // Test C — Quiescent PSC produces nothing.
    Cell* psc_c = create_cell(*pStroma);
    psc_c->assign_position(std::vector<double>{140.0, 100.0, 0.0});
    const int c_acta2 = psc_c->custom_data.find_variable_index("acta2_active");
    assert(c_acta2 >= 0);
    psc_c->custom_data[c_acta2] = 0.0;
    const int voxel_c = voxel_index_for_cell(psc_c);
    set_ecm_at_voxel(voxel_c, 0.3, 0.5);
    module6_ecm_production(psc_c, psc_c->phenotype, dt, ModulePhase::WRITE);
    assert(nearly_equal(get_ecm_density_at_voxel(voxel_c), 0.3));
    std::cout << "PASS Test C" << std::endl;

    // Test D — MMP2 degrades ECM proportionally.
    Cell* tumor_d = create_cell(*pTumor);
    tumor_d->assign_position(std::vector<double>{180.0, 100.0, 0.0});
    const int d_mmp2 = tumor_d->custom_data.find_variable_index("mmp2_active");
    assert(d_mmp2 >= 0);
    tumor_d->custom_data[d_mmp2] = 1.0;
    parameters.doubles("mmp2_degradation_rate") = 0.005;
    const int voxel_d = voxel_index_for_cell(tumor_d);
    set_ecm_at_voxel(voxel_d, 0.5, 0.6);
    module6_ecm_production(tumor_d, tumor_d->phenotype, dt, ModulePhase::WRITE);
    assert(nearly_equal(get_ecm_density_at_voxel(voxel_d), 0.47));
    assert(nearly_equal(get_ecm_ha_fraction(voxel_d), 0.6));
    std::cout << "PASS Test D" << std::endl;

    // Test E — ECM density never goes below 0.
    parameters.doubles("mmp2_degradation_rate") = 0.1;
    set_ecm_at_voxel(voxel_d, 0.01, 0.6);
    module6_ecm_production(tumor_d, tumor_d->phenotype, dt, ModulePhase::WRITE);
    assert(nearly_equal(get_ecm_density_at_voxel(voxel_d), 0.0));
    std::cout << "PASS Test E" << std::endl;

    // Test F — ECM density never goes above 1.
    Cell* caf_f = create_cell(*pStroma);
    caf_f->assign_position(std::vector<double>{220.0, 100.0, 0.0});
    const int f_acta2 = caf_f->custom_data.find_variable_index("acta2_active");
    const int f_gli1 = caf_f->custom_data.find_variable_index("gli1_active");
    assert(f_acta2 >= 0);
    assert(f_gli1 >= 0);
    caf_f->custom_data[f_acta2] = 1.0;
    caf_f->custom_data[f_gli1] = 1.0;
    parameters.doubles("ecm_production_rate_boosted") = 0.1;
    const int voxel_f = voxel_index_for_cell(caf_f);
    set_ecm_at_voxel(voxel_f, 0.99, 0.5);
    module6_ecm_production(caf_f, caf_f->phenotype, dt, ModulePhase::WRITE);
    assert(nearly_equal(get_ecm_density_at_voxel(voxel_f), 1.0));
    std::cout << "PASS Test F" << std::endl;

    // Test G — HA degrade only (independence).
    Cell* cell_g = create_cell(*pStroma);
    cell_g->assign_position(std::vector<double>{260.0, 100.0, 0.0});
    const int g_acta2 = cell_g->custom_data.find_variable_index("acta2_active");
    assert(g_acta2 >= 0);
    cell_g->custom_data[g_acta2] = 0.0; // no production
    const int voxel_g = voxel_index_for_cell(cell_g);
    set_ecm_at_voxel(voxel_g, 0.5, 0.6); // HA=0.3, COL=0.2
    intervention_state.ha_degrade_active = true;
    intervention_state.col_degrade_active = false;
    intervention_state.ha_degrade_strength = 0.01;
    intervention_state.col_degrade_strength = 0.0;
    module6_ecm_production(cell_g, cell_g->phenotype, dt, ModulePhase::WRITE);
    const double density_g = get_ecm_density_at_voxel(voxel_g);
    const double frac_g = get_ecm_ha_fraction(voxel_g);
    const double col_g = density_g * (1.0 - frac_g);
    assert(nearly_equal(density_g, 0.44));
    assert(nearly_equal(col_g, 0.2));
    std::cout << "PASS Test G" << std::endl;

    // Test H — COL degrade only (independence).
    Cell* cell_h = create_cell(*pStroma);
    cell_h->assign_position(std::vector<double>{300.0, 100.0, 0.0});
    const int h_acta2 = cell_h->custom_data.find_variable_index("acta2_active");
    assert(h_acta2 >= 0);
    cell_h->custom_data[h_acta2] = 0.0; // no production
    const int voxel_h = voxel_index_for_cell(cell_h);
    set_ecm_at_voxel(voxel_h, 0.5, 0.6); // HA=0.3, COL=0.2
    intervention_state.ha_degrade_active = false;
    intervention_state.col_degrade_active = true;
    intervention_state.ha_degrade_strength = 0.0;
    intervention_state.col_degrade_strength = 0.01;
    module6_ecm_production(cell_h, cell_h->phenotype, dt, ModulePhase::WRITE);
    const double density_h = get_ecm_density_at_voxel(voxel_h);
    const double frac_h = get_ecm_ha_fraction(voxel_h);
    const double ha_h = density_h * frac_h;
    assert(nearly_equal(density_h, 0.44));
    assert(nearly_equal(ha_h, 0.3));
    std::cout << "PASS Test H" << std::endl;

    intervention_state.ha_degrade_active = false;
    intervention_state.col_degrade_active = false;
    intervention_state.ha_degrade_strength = 0.0;
    intervention_state.col_degrade_strength = 0.0;

    std::cout << "PASS module6_ecm_production_test" << std::endl;
    return 0;
}
