#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "feedback_loop_common.h"

using namespace BioFVM;
using namespace PhysiCell;
using namespace feedback_loop_common;

namespace
{

double mean_field_in_annulus_local(int substrate_index,
                                   double cx,
                                   double cy,
                                   double r_inner,
                                   double r_outer)
{
    const int n_vox = static_cast<int>(microenvironment.number_of_voxels());
    double sum = 0.0;
    int count = 0;
    for (int n = 0; n < n_vox; ++n)
    {
        const std::vector<double>& c = microenvironment.mesh.voxels[n].center;
        const double dx = c[0] - cx;
        const double dy = c[1] - cy;
        const double r = std::sqrt(dx * dx + dy * dy);
        if (r >= r_inner && r <= r_outer)
        {
            sum += microenvironment.density_vector(n)[substrate_index];
            ++count;
        }
    }
    return (count > 0) ? (sum / static_cast<double>(count)) : 0.0;
}

void set_uniform_substrate(int substrate_index, double value)
{
    const int n_vox = static_cast<int>(microenvironment.number_of_voxels());
    for (int n = 0; n < n_vox; ++n)
    {
        microenvironment.density_vector(n)[substrate_index] = value;
    }
}

void set_uniform_ecm(double density, double ha_fraction)
{
    const int n_vox = static_cast<int>(microenvironment.number_of_voxels());
    for (int n = 0; n < n_vox; ++n)
    {
        microenvironment.density_vector(n)[ecm_index] = density;
        set_ecm_ha_fraction(n, ha_fraction);
    }
}

struct Metrics
{
    double at_source = 0.0;
    double beyond_barrier = 0.0;
};

Metrics run_condition(bool with_ecm_ring)
{
    set_uniform_substrate(tgfb_index, 0.0);
    set_uniform_ecm(0.0, 0.6);

    if (with_ecm_ring)
    {
        set_annulus_ecm(0.0, 0.0, 60.0, 180.0, 0.8, 0.8);
    }

    const int center_voxel = voxel_index_for_position(0.0, 0.0, 0.0);

    const double dt = 1.0;
    for (int step = 0; step < 200; ++step)
    {
        // Continuous local source to mimic persistent tumor secretion.
        microenvironment.density_vector(center_voxel)[tgfb_index] = 1.0;
        microenvironment.simulate_diffusion_decay(dt);
    }

    Metrics m;
    m.at_source = mean_field_in_annulus_local(tgfb_index, 0.0, 0.0, 0.0, 20.0);
    m.beyond_barrier = mean_field_in_annulus_local(tgfb_index, 0.0, 0.0, 260.0, 380.0);
    return m;
}

} // namespace

int main()
{
    initialize_world();

    assert(tgfb_index >= 0);
    assert(ecm_index >= 0);

    const Metrics no_ecm = run_condition(false);
    const Metrics with_ecm = run_condition(true);

    std::cout << "Rule36 metrics"
              << " source_no_ecm=" << no_ecm.at_source
              << " source_with_ecm=" << with_ecm.at_source
              << " beyond_no_ecm=" << no_ecm.beyond_barrier
              << " beyond_with_ecm=" << with_ecm.beyond_barrier
              << std::endl;

    // Robust trapping signature in this solver configuration: less escaped far away.
    assert(with_ecm.beyond_barrier < no_ecm.beyond_barrier);

    std::cout << "PASS Rule36_ecm_traps_paracrine_signals" << std::endl;
    std::cout << "PASS rule36_ecm_traps_paracrine_test" << std::endl;
    return 0;
}
