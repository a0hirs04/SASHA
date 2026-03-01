#include <cassert>
#include <cmath>
#include <iostream>

#include "../Stroma_world/PhysiCell/core/PhysiCell.h"
#include "../custom_modules/custom.h"

using namespace BioFVM;
using namespace PhysiCell;

namespace
{

bool nearly_equal(double a, double b, double eps = 1e-10)
{
    return std::fabs(a - b) <= eps;
}

void set_voxel_ecm_state(int voxel_index, double density, double ha_fraction)
{
    std::vector<double>& rho = microenvironment.density_vector(voxel_index);
    rho[ecm_index] = density;
    set_ecm_ha_fraction(voxel_index, ha_fraction);
    update_ecm_effective_diffusion_coefficients(microenvironment);
}

} // namespace

int main()
{
    const bool xml_status = load_PhysiCell_config_file("config/PhysiCell_settings.xml");
    assert(xml_status);

    setup_microenvironment();

    assert(oxygen_index >= 0);
    assert(tgfb_index >= 0);
    assert(shh_index >= 0);
    assert(drug_index >= 0);
    assert(ecm_index >= 0);

    std::vector<double> probe_position{100.0, 100.0, 0.0};
    const int voxel = microenvironment.nearest_voxel_index(probe_position);
    assert(voxel >= 0);

    const double d_o2_base = microenvironment.diffusion_coefficients[oxygen_index];
    const double d_tgfb_base = microenvironment.diffusion_coefficients[tgfb_index];
    const double d_shh_base = microenvironment.diffusion_coefficients[shh_index];
    const double d_drug_base = microenvironment.diffusion_coefficients[drug_index];

    // Test A — No ECM, full diffusion.
    set_voxel_ecm_state(voxel, 0.0, 0.5);
    assert(nearly_equal(get_effective_diffusion_coefficient(oxygen_index, voxel), d_o2_base));
    assert(nearly_equal(get_effective_diffusion_coefficient(tgfb_index, voxel), d_tgfb_base));
    assert(nearly_equal(get_effective_diffusion_coefficient(shh_index, voxel), d_shh_base));
    assert(nearly_equal(get_effective_diffusion_coefficient(drug_index, voxel), d_drug_base));
    std::cout << "PASS Test A" << std::endl;

    // Rule 33 — ECM monotonically decreases drug diffusion.
    set_voxel_ecm_state(voxel, 0.0, 0.5);
    const double d0 = get_effective_diffusion_coefficient(drug_index, voxel);
    set_voxel_ecm_state(voxel, 0.3, 0.5);
    const double d03 = get_effective_diffusion_coefficient(drug_index, voxel);
    set_voxel_ecm_state(voxel, 0.6, 0.5);
    const double d06 = get_effective_diffusion_coefficient(drug_index, voxel);
    set_voxel_ecm_state(voxel, 0.9, 0.5);
    const double d09 = get_effective_diffusion_coefficient(drug_index, voxel);
    assert(d0 > d03 && d03 > d06 && d06 > d09);
    std::cout << "PASS Rule33_monotonic_diffusion_impedance" << std::endl;

    // Test B — Dense ECM, HA-dominant.
    set_voxel_ecm_state(voxel, 0.8, 0.8);
    const double d_drug_b = get_effective_diffusion_coefficient(drug_index, voxel);
    const double d_drug_b_expected = d_drug_base * 0.44; // 1 - 0.8 * (0.8*0.8 + 0.3*0.2)
    assert(nearly_equal(d_drug_b, d_drug_b_expected));
    std::cout << "PASS Test B" << std::endl;

    // Test C — Same density, collagen-dominant.
    set_voxel_ecm_state(voxel, 0.8, 0.2);
    const double d_drug_c = get_effective_diffusion_coefficient(drug_index, voxel);
    const double d_drug_c_expected = d_drug_base * 0.68; // 1 - 0.8 * (0.8*0.2 + 0.3*0.8)
    assert(nearly_equal(d_drug_c, d_drug_c_expected));
    assert(d_drug_c > d_drug_b);
    std::cout << "PASS Test C" << std::endl;
    std::cout << "PASS Rule34_ha_impedes_more_than_collagen" << std::endl;

    // Test D — HA depletion improves drug diffusion.
    set_voxel_ecm_state(voxel, 0.8, 0.8);
    const double d_drug_before_ha_depletion = get_effective_diffusion_coefficient(drug_index, voxel);
    set_voxel_ecm_state(voxel, 0.8, 0.2);
    const double d_drug_after_ha_depletion = get_effective_diffusion_coefficient(drug_index, voxel);
    assert(nearly_equal(d_drug_after_ha_depletion, d_drug_c_expected));
    assert(d_drug_after_ha_depletion > d_drug_before_ha_depletion);
    std::cout << "PASS Test D" << std::endl;

    // Test E — Collagen depletion changes drug diffusion less than HA depletion.
    set_voxel_ecm_state(voxel, 0.8, 0.8);
    const double d_drug_reference = get_effective_diffusion_coefficient(drug_index, voxel);
    set_voxel_ecm_state(voxel, 0.8, 0.9); // collagen depletion proxy
    const double d_drug_after_col_depletion = get_effective_diffusion_coefficient(drug_index, voxel);
    set_voxel_ecm_state(voxel, 0.8, 0.2); // HA depletion
    const double d_drug_after_ha_depletion_e = get_effective_diffusion_coefficient(drug_index, voxel);
    const double col_change = std::fabs(d_drug_after_col_depletion - d_drug_reference);
    const double ha_change = std::fabs(d_drug_after_ha_depletion_e - d_drug_reference);
    assert(col_change < ha_change);
    std::cout << "PASS Test E" << std::endl;

    // Rule 35 — collagen-dominant ECM is mechanically stiffer than HA-dominant ECM.
    set_voxel_ecm_state(voxel, 0.7, 0.8); // HA-dominant (collagen 0.2)
    const double stiff_ha_dom = get_local_mechanical_stiffness(voxel);
    set_voxel_ecm_state(voxel, 0.7, 0.2); // collagen-dominant (collagen 0.8)
    const double stiff_col_dom = get_local_mechanical_stiffness(voxel);
    assert(stiff_col_dom > stiff_ha_dom);
    std::cout << "PASS Rule35_collagen_stiffer_than_HA" << std::endl;

    std::cout << "PASS ecm_diffusion_coupling_test" << std::endl;
    return 0;
}
