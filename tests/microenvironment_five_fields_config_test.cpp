#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../Stroma_world/PhysiCell/core/PhysiCell.h"
#include "../custom_modules/custom.h"

using namespace BioFVM;
using namespace PhysiCell;

namespace
{

bool nearly_equal(double a, double b, double tol = 1e-8)
{
    return std::fabs(a - b) <= tol;
}

bool is_edge_voxel(const std::vector<double>& c,
                   double xmin, double xmax,
                   double ymin, double ymax,
                   double dx, double dy)
{
    const double xlo = xmin + 0.5 * dx;
    const double xhi = xmax - 0.5 * dx;
    const double ylo = ymin + 0.5 * dy;
    const double yhi = ymax - 0.5 * dy;
    return (c[0] <= xlo + 1e-12) || (c[0] >= xhi - 1e-12) ||
           (c[1] <= ylo + 1e-12) || (c[1] >= yhi - 1e-12);
}

} // namespace

int main()
{
    const bool ok = load_PhysiCell_config_file("config/PhysiCell_settings.xml");
    assert(ok);
    setup_microenvironment();

    // 1) Exactly 5 fields.
    assert(microenvironment.density_names.size() == 5);
    assert(oxygen_index >= 0 && tgfb_index >= 0 && shh_index >= 0 && drug_index >= 0 && ecm_index >= 0);
    assert(microenvironment.find_density_index("oxygen") == oxygen_index);
    assert(microenvironment.find_density_index("tgfb") == tgfb_index);
    assert(microenvironment.find_density_index("shh") == shh_index);
    assert(microenvironment.find_density_index("drug") == drug_index);
    assert(microenvironment.find_density_index("ecm_density") == ecm_index);

    // 2) ECM does not diffuse.
    assert(nearly_equal(microenvironment.diffusion_coefficients[ecm_index], 0.0));

    const double xmin = microenvironment.mesh.bounding_box[0];
    const double ymin = microenvironment.mesh.bounding_box[1];
    const double xmax = microenvironment.mesh.bounding_box[3];
    const double ymax = microenvironment.mesh.bounding_box[4];
    const double dx = microenvironment.mesh.dx;
    const double dy = microenvironment.mesh.dy;

    int edge_count = 0;
    for (unsigned int n = 0; n < microenvironment.number_of_voxels(); ++n)
    {
        const std::vector<double>& c = microenvironment.mesh.voxels[n].center;
        if (!is_edge_voxel(c, xmin, xmax, ymin, ymax, dx, dy)) continue;
        ++edge_count;

        const std::vector<double>& rho = microenvironment.density_vector(static_cast<int>(n));
        // 3) O2 boundary fixed.
        assert(nearly_equal(rho[oxygen_index], 38.0, 1e-6));
        // 4) Drug boundary baseline 0.
        assert(nearly_equal(rho[drug_index], 0.0, 1e-12));
        // 5) TGFb / SHH boundary sinks at 0.
        assert(nearly_equal(rho[tgfb_index], 0.0, 1e-12));
        assert(nearly_equal(rho[shh_index], 0.0, 1e-12));
        // 6) ECM boundary 0.
        assert(nearly_equal(rho[ecm_index], 0.0, 1e-12));
    }
    assert(edge_count > 0);

    // Choose an interior voxel.
    int interior_voxel = -1;
    for (unsigned int n = 0; n < microenvironment.number_of_voxels(); ++n)
    {
        const std::vector<double>& c = microenvironment.mesh.voxels[n].center;
        if (!is_edge_voxel(c, xmin, xmax, ymin, ymax, dx, dy))
        {
            interior_voxel = static_cast<int>(n);
            break;
        }
    }
    assert(interior_voxel >= 0);

    // 8) ecm_ha_fraction custom attribute is per-voxel and accessible.
    set_ecm_ha_fraction(interior_voxel, 0.6);
    assert(nearly_equal(get_ecm_ha_fraction(interior_voxel), 0.6, 1e-12));
    assert(nearly_equal(get_ecm_collagen_fraction(interior_voxel), 0.4, 1e-12));

    // 7) ECM remains local (no diffusion), only tiny decay is allowed.
    const std::vector<double>& c0 = microenvironment.mesh.voxels[interior_voxel].center;
    std::vector<double> c1 = c0;
    c1[0] += dx; // adjacent in x
    int adjacent_voxel = microenvironment.nearest_voxel_index(c1);
    if (adjacent_voxel == interior_voxel)
    {
        c1 = c0;
        c1[1] += dy; // fallback adjacent in y
        adjacent_voxel = microenvironment.nearest_voxel_index(c1);
    }
    assert(adjacent_voxel >= 0);
    assert(adjacent_voxel != interior_voxel);

    microenvironment.density_vector(interior_voxel)[ecm_index] = 0.8;
    microenvironment.density_vector(adjacent_voxel)[ecm_index] = 0.0;

    const double dt = 1.0; // min
    for (int i = 0; i < 100; ++i)
    {
        microenvironment.simulate_diffusion_decay(dt);
    }
    const double ecm_after_interior = microenvironment.density_vector(interior_voxel)[ecm_index];
    const double ecm_after_adjacent = microenvironment.density_vector(adjacent_voxel)[ecm_index];
    assert(std::fabs(ecm_after_interior - 0.8) < 0.01); // approximately unchanged (slow decay only)
    assert(std::fabs(ecm_after_adjacent - 0.0) < 1e-10); // no spill into adjacent voxel
    std::cout << "PASS Rule28_ecm_non_diffusive" << std::endl;

    // 29) ECM near-permanent: 500 steps, <5% decrease with no producers/degraders.
    microenvironment.density_vector(interior_voxel)[ecm_index] = 0.5;
    const double ecm_start_29 = microenvironment.density_vector(interior_voxel)[ecm_index];
    for (int i = 0; i < 500; ++i)
    {
        microenvironment.simulate_diffusion_decay(dt);
    }
    const double ecm_end_29 = microenvironment.density_vector(interior_voxel)[ecm_index];
    const double frac_drop_29 = (ecm_start_29 - ecm_end_29) / ecm_start_29;
    assert(frac_drop_29 < 0.05);
    std::cout << "PASS Rule29_ecm_near_permanent" << std::endl;

    // Derived properties should be available after diffusion update.
    const double o2_eff = get_effective_diffusion_coefficient(oxygen_index, interior_voxel);
    const double mech_stiff = get_local_mechanical_stiffness(interior_voxel);
    assert(o2_eff > 0.0);
    assert(mech_stiff >= 0.0);

    // 9) Oxygen nowhere exceeds boundary value.
    double max_o2 = -1e99;
    for (unsigned int n = 0; n < microenvironment.number_of_voxels(); ++n)
    {
        const std::vector<double>& rho = microenvironment.density_vector(static_cast<int>(n));
        if (rho[oxygen_index] > max_o2) max_o2 = rho[oxygen_index];
        // Rule 63 clamp range for ECM
        assert(rho[ecm_index] >= -1e-12);
        assert(rho[ecm_index] <= 1.0 + 1e-12);
    }
    assert(max_o2 <= 38.0 + 1e-8);

    std::cout << "PASS microenvironment_five_fields_config_test" << std::endl;
    return 0;
}
