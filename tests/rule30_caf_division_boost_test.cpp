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

int count_stromal_window(int stromal_type, double x_min, double x_max)
{
    int n = 0;
    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->phenotype.death.dead) continue;
        if (pCell->type != stromal_type) continue;
        if (pCell->position[0] >= x_min && pCell->position[0] <= x_max) ++n;
    }
    return n;
}

void set_local_signal(Cell* pCell, int substrate_index, double value)
{
    std::vector<double>& rho = pCell->nearest_density_vector();
    if (substrate_index >= 0 && substrate_index < static_cast<int>(rho.size()))
    {
        rho[substrate_index] = value;
    }
}

} // namespace

int main()
{
    initialize_world();

    Cell_Definition* pStroma = find_cell_definition("stromal_cell");
    assert(pStroma != NULL);
    assert(tgfb_index >= 0);
    assert(shh_index >= 0);

    // Force moderate division with clear GLI1 boost signal.
    parameters.doubles("caf_proliferation_rate") = 0.001;
    parameters.doubles("gli1_proliferation_boost") = 2.0;
    parameters.doubles("psc_proliferation_rate") = 0.0;
    parameters.doubles("contact_inhibition_threshold") = 1e9;
    parameters.doubles("stromal_contact_inhibition_threshold") = 1e9;
    parameters.doubles("caf_base_death_rate") = 0.0;
    parameters.doubles("caf_crowded_death_rate") = 0.0;
    parameters.doubles("ecm_production_rate_base") = 0.0;
    parameters.doubles("ecm_production_rate_boosted") = 0.0;

    // Keep fields simple.
    reset_all_fields(38.0, 0.0, 0.0, 0.0, 0.0, 0.5);

    Cell* rep_on = NULL;
    Cell* rep_off = NULL;

    // Left: GLI1-ON CAFs. Right: GLI1-OFF CAFs.
    for (int i = 0; i < 40; ++i)
    {
        Cell* on = create_cell(*pStroma);
        on->assign_position(std::vector<double>{-420.0 + 8.0 * i, 0.0, 0.0});
        on->custom_data[custom_index(on, "acta2_active")] = 1.0;
        on->custom_data[custom_index(on, "activation_mode")] = 1.0;
        on->custom_data[custom_index(on, "gli1_active")] = 1.0;
        on->phenotype.motility.is_motile = false;

        // Maintain GLI1 ON via SHH.
        set_local_signal(on, tgfb_index, 0.0);
        set_local_signal(on, shh_index, 0.5);
        if (rep_on == NULL) rep_on = on;

        Cell* off = create_cell(*pStroma);
        off->assign_position(std::vector<double>{120.0 + 8.0 * i, 0.0, 0.0});
        off->custom_data[custom_index(off, "acta2_active")] = 1.0;
        off->custom_data[custom_index(off, "activation_mode")] = 1.0;
        off->custom_data[custom_index(off, "gli1_active")] = 0.0;
        off->phenotype.motility.is_motile = false;

        // Maintain GLI1 OFF (no SHH).
        set_local_signal(off, tgfb_index, 0.0);
        set_local_signal(off, shh_index, 0.0);
        if (rep_off == NULL) rep_off = off;
    }

    const int stromal_type = pStroma->type;
    const double on_x_min = -520.0;
    const double on_x_max = -40.0;
    const double off_x_min = 40.0;
    const double off_x_max = 520.0;
    const int initial_on = count_stromal_window(stromal_type, on_x_min, on_x_max);
    const int initial_off = count_stromal_window(stromal_type, off_x_min, off_x_max);
    assert(initial_on == 40);
    assert(initial_off == 40);

    double t = 0.0;
    const double dt = 6.0;
    Cell_Container* cc = static_cast<Cell_Container*>(microenvironment.agent_container);
    for (int step = 0; step < 80; ++step)
    {
        for (size_t i = 0; i < all_cells->size(); ++i)
        {
            Cell* pCell = (*all_cells)[i];
            if (pCell == NULL || pCell->phenotype.death.dead) continue;
            if (pCell->type != stromal_type) continue;

            // Keep SHH split so GLI1 remains ON only in the left cohort.
            if (pCell->position[0] < 0.0)
            {
                set_local_signal(pCell, shh_index, 0.5);
                pCell->custom_data[custom_index(pCell, "gli1_active")] = 1.0;
            }
            else
            {
                set_local_signal(pCell, shh_index, 0.0);
                pCell->custom_data[custom_index(pCell, "gli1_active")] = 0.0;
            }
            set_local_signal(pCell, tgfb_index, 0.0);
        }

        microenvironment.simulate_diffusion_decay(dt);
        cc->update_all_cells(t, dt, dt, dt);
        t += dt;
    }

    const int final_on = count_stromal_window(stromal_type, on_x_min, on_x_max);
    const int final_off = count_stromal_window(stromal_type, off_x_min, off_x_max);

    assert(rep_on != NULL);
    assert(rep_off != NULL);
    const double rate_on = rep_on->phenotype.cycle.data.transition_rate(0, 0);
    const double rate_off = rep_off->phenotype.cycle.data.transition_rate(0, 0);

    std::cout << "Rule30 metrics"
              << " initial_on=" << initial_on
              << " initial_off=" << initial_off
              << " final_on=" << final_on
              << " final_off=" << final_off
              << " rate_on=" << rate_on
              << " rate_off=" << rate_off
              << std::endl;

    assert(final_on > initial_on);
    assert(final_off > initial_off);
    assert(rate_on > rate_off);

    std::cout << "PASS Rule30_caf_division_gli1_boost" << std::endl;
    std::cout << "PASS rule30_caf_division_boost_test" << std::endl;
    return 0;
}
