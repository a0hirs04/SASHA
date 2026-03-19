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

bool nearly_equal(double a, double b, double eps = 1e-6)
{
    return std::fabs(a - b) <= eps;
}

void set_local_signal(Cell* pCell, int substrate_index, double value)
{
    std::vector<double>& rho = pCell->nearest_density_vector();
    if (substrate_index >= 0 && substrate_index < static_cast<int>(rho.size()))
    {
        rho[substrate_index] = value;
    }
}

int count_live_stromal_in_x_window(int stromal_type, double x_min, double x_max)
{
    int n = 0;
    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->phenotype.death.dead) continue;
        if (pCell->type != stromal_type) continue;
        const double x = pCell->position[0];
        if (x >= x_min && x <= x_max) ++n;
    }
    return n;
}

} // namespace

int main()
{
    initialize_world();

    Cell_Definition* pStroma = find_cell_definition("stromal_cell");
    assert(pStroma != NULL);
    assert(ecm_index >= 0);

    // Shared isolation.
    parameters.doubles("base_proliferation_rate") = 0.0;
    parameters.doubles("apoptosis_resistance") = 1.0;
    parameters.doubles("drug_kill_coefficient") = 0.0;
    parameters.doubles("efflux_drug_reduction") = 0.0;
    parameters.doubles("tgfb_secretion_rate") = 0.0;
    parameters.doubles("shh_secretion_rate") = 0.0;
    parameters.doubles("mechanical_compaction_strength") = 0.0;
    parameters.doubles("compaction_ecm_increment") = 0.0;
    parameters.doubles("mmp2_degradation_rate") = 0.0;

    reset_all_fields(0.25, 0.0, 0.0, 0.0, 0.0, 0.6);

    // E11: ECM production gli1 ON > OFF in activated CAFs.
    parameters.doubles("ecm_production_rate_base") = 0.01;
    parameters.doubles("ecm_production_rate_boosted") = 0.015;
    parameters.doubles("ecm_ha_fraction_default") = 0.6;

    Cell* caf_off = create_cell(*pStroma);
    caf_off->assign_position(std::vector<double>{-150.0, 0.0, 0.0});
    caf_off->phenotype.motility.is_motile = false;
    Cell* caf_on = create_cell(*pStroma);
    caf_on->assign_position(std::vector<double>{150.0, 0.0, 0.0});
    caf_on->phenotype.motility.is_motile = false;

    caf_off->custom_data[custom_index(caf_off, "acta2_active")] = 1.0;
    caf_off->custom_data[custom_index(caf_off, "gli1_active")] = 0.0;
    caf_on->custom_data[custom_index(caf_on, "acta2_active")] = 1.0;
    caf_on->custom_data[custom_index(caf_on, "gli1_active")] = 1.0;
    set_local_signal(caf_off, tgfb_index, 0.0);
    set_local_signal(caf_off, shh_index, 0.0);
    set_local_signal(caf_on, tgfb_index, 0.0);
    set_local_signal(caf_on, shh_index, 0.5);

    const int voxel_off = voxel_index_for_cell(caf_off);
    const int voxel_on = voxel_index_for_cell(caf_on);
    microenvironment.density_vector(voxel_off)[ecm_index] = 0.0;
    microenvironment.density_vector(voxel_on)[ecm_index] = 0.0;

    double t = 0.0;
    const double dt = 1.0;
    for (int step = 0; step < 50; ++step)
    {
        set_local_signal(caf_off, tgfb_index, 0.0);
        set_local_signal(caf_off, shh_index, 0.0);
        set_local_signal(caf_on, tgfb_index, 0.0);
        set_local_signal(caf_on, shh_index, 0.5);
        advance_steps(1, dt, t);
    }

    const double ecm_off = microenvironment.density_vector(voxel_off)[ecm_index];
    const double ecm_on = microenvironment.density_vector(voxel_on)[ecm_index];
    const bool pass_e11 = (ecm_on > ecm_off);

    // E13: CAF proliferation boost via GLI1.
    parameters.doubles("caf_proliferation_rate") = 0.003;
    parameters.doubles("gli1_proliferation_boost") = 2.0;
    parameters.doubles("psc_proliferation_rate") = 0.0;

    Cell* rep_off = NULL;
    Cell* rep_on = NULL;

    const double two_pi = 6.283185307179586;
    for (int i = 0; i < 5; ++i)
    {
        const double th = two_pi * (static_cast<double>(i) / 5.0);
        Cell* c0 = create_cell(*pStroma);
        c0->assign_position(std::vector<double>{-500.0 + 20.0 * std::cos(th), 120.0 + 20.0 * std::sin(th), 0.0});
        c0->phenotype.motility.is_motile = false;
        c0->custom_data[custom_index(c0, "acta2_active")] = 1.0;
        c0->custom_data[custom_index(c0, "gli1_active")] = 0.0;
        set_local_signal(c0, tgfb_index, 0.0);
        set_local_signal(c0, shh_index, 0.0);
        if (rep_off == NULL) rep_off = c0;

        Cell* c1 = create_cell(*pStroma);
        c1->assign_position(std::vector<double>{500.0 + 20.0 * std::cos(th), 120.0 + 20.0 * std::sin(th), 0.0});
        c1->phenotype.motility.is_motile = false;
        c1->custom_data[custom_index(c1, "acta2_active")] = 1.0;
        c1->custom_data[custom_index(c1, "gli1_active")] = 1.0;
        set_local_signal(c1, tgfb_index, 0.0);
        set_local_signal(c1, shh_index, 0.5);
        if (rep_on == NULL) rep_on = c1;
    }

    // Count cohorts by spatial windows to include daughter cells.
    const double off_x_min = -900.0;
    const double off_x_max = -300.0;
    const double on_x_min = 300.0;
    const double on_x_max = 900.0;

    const int init_off = count_live_stromal_in_x_window(pStroma->type, off_x_min, off_x_max);
    const int init_on = count_live_stromal_in_x_window(pStroma->type, on_x_min, on_x_max);

    for (int step = 0; step < 1; ++step)
    {
        for (Cell* pCell : *all_cells)
        {
            if (pCell == NULL || pCell->phenotype.death.dead) continue;
            if (pCell->type != pStroma->type) continue;
            if (pCell->position[0] >= on_x_min && pCell->position[0] <= on_x_max)
            {
                set_local_signal(pCell, tgfb_index, 0.0);
                set_local_signal(pCell, shh_index, 0.5);
            }
            else if (pCell->position[0] >= off_x_min && pCell->position[0] <= off_x_max)
            {
                set_local_signal(pCell, tgfb_index, 0.0);
                set_local_signal(pCell, shh_index, 0.0);
            }
        }
        advance_steps(1, dt, t);
    }

    assert(rep_off != NULL);
    assert(rep_on != NULL);

    // Requested diagnostics: cycle model and transition indices/rates at runtime.
    const Cycle_Model& stromal_cycle_model = rep_off->phenotype.cycle.model();
    const std::string cycle_model_name = stromal_cycle_model.name;
    const int configured_start_index = 0;
    const int configured_end_index = 0; // module4 writes transition_rate(0,0)

    int runtime_start_index = 0;
    int runtime_end_index = -1;
    if (!stromal_cycle_model.phase_links.empty() &&
        !stromal_cycle_model.phase_links[0].empty())
    {
        runtime_end_index = stromal_cycle_model.phase_links[0][0].end_phase_index;
    }

    const double configured_rate_off =
        rep_off->phenotype.cycle.data.transition_rate(configured_start_index, configured_end_index);
    const double configured_rate_on =
        rep_on->phenotype.cycle.data.transition_rate(configured_start_index, configured_end_index);

    double runtime_rate_off = configured_rate_off;
    double runtime_rate_on = configured_rate_on;
    if (runtime_end_index >= 0)
    {
        runtime_rate_off =
            rep_off->phenotype.cycle.data.transition_rate(runtime_start_index, runtime_end_index);
        runtime_rate_on =
            rep_on->phenotype.cycle.data.transition_rate(runtime_start_index, runtime_end_index);
    }

    for (int step = 0; step < 199; ++step)
    {
        for (Cell* pCell : *all_cells)
        {
            if (pCell == NULL || pCell->phenotype.death.dead) continue;
            if (pCell->type != pStroma->type) continue;
            if (pCell->position[0] >= on_x_min && pCell->position[0] <= on_x_max)
            {
                set_local_signal(pCell, tgfb_index, 0.0);
                set_local_signal(pCell, shh_index, 0.5);
            }
            else if (pCell->position[0] >= off_x_min && pCell->position[0] <= off_x_max)
            {
                set_local_signal(pCell, tgfb_index, 0.0);
                set_local_signal(pCell, shh_index, 0.0);
            }
        }
        advance_steps(1, dt, t);
    }

    const int final_off = count_live_stromal_in_x_window(pStroma->type, off_x_min, off_x_max);
    const int final_on = count_live_stromal_in_x_window(pStroma->type, on_x_min, on_x_max);

    const double expected_boost = parameters.doubles("gli1_proliferation_boost");
    const double rate_ratio = (runtime_rate_off > 0.0) ? (runtime_rate_on / runtime_rate_off) : 0.0;

    const bool pass_e13_ratio = nearly_equal(rate_ratio, expected_boost);
    // E13 passes only on actual population divergence.
    const bool pass_e13_count = (final_on > final_off);

    const bool pass = pass_e11 && pass_e13_count;

    std::cout << "E13 diagnostics"
              << " cycle_model_name=\"" << cycle_model_name << "\""
              << " configured_transition_index=" << configured_start_index << "->" << configured_end_index
              << " runtime_transition_index=" << runtime_start_index << "->" << runtime_end_index
              << " configured_rate_off=" << configured_rate_off
              << " configured_rate_on=" << configured_rate_on
              << " runtime_rate_off=" << runtime_rate_off
              << " runtime_rate_on=" << runtime_rate_on
              << std::endl;

    std::cout << "P5B measurements"
              << " ecm_gli1_off=" << ecm_off
              << " ecm_gli1_on=" << ecm_on
              << " runtime_rate_off=" << runtime_rate_off
              << " runtime_rate_on=" << runtime_rate_on
              << " rate_ratio=" << rate_ratio
              << " expected_boost=" << expected_boost
              << " init_off_n=" << init_off
              << " init_on_n=" << init_on
              << " final_off_n=" << final_off
              << " final_on_n=" << final_on
              << std::endl;

    std::cout << "P5B checks"
              << " E11_ecm_on_gt_off=" << (pass_e11 ? 1 : 0)
              << " E13_rate_ratio_matches_boost=" << (pass_e13_ratio ? 1 : 0)
              << " E13_count_on_gt_off=" << (pass_e13_count ? 1 : 0)
              << std::endl;

    if (!pass)
    {
        std::cout << "FAIL P5 Cluster B" << std::endl;
        return 1;
    }

    std::cout << "PASS P5 Cluster B" << std::endl;
    return 0;
}
