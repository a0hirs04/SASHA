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

double max_field_value(int substrate_index)
{
    double mx = 0.0;
    const int n_vox = static_cast<int>(microenvironment.number_of_voxels());
    for (int n = 0; n < n_vox; ++n)
    {
        mx = std::max(mx, microenvironment.density_vector(n)[substrate_index]);
    }
    return mx;
}

} // namespace

int main()
{
    initialize_world();

    assert(tgfb_index >= 0);
    assert(shh_index >= 0);

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    Cell_Definition* pStroma = find_cell_definition("stromal_cell");
    assert(pTumor != NULL);
    assert(pStroma != NULL);

    // E07 setup: SHH-only arm.
    parameters.doubles("tgfb_secretion_rate") = 0.0;
    parameters.doubles("shh_secretion_rate") = 0.5;
    parameters.doubles("hif1a_tgfb_amplification_factor") = 1.0;
    parameters.doubles("caf_tgfb_secretion_rate") = 0.0;

    parameters.doubles("tgfb_activation_threshold") = 0.3;
    parameters.doubles("shh_activation_threshold") = 0.05;

    parameters.doubles("base_proliferation_rate") = 0.0;
    parameters.doubles("caf_proliferation_rate") = 0.0;
    parameters.doubles("psc_proliferation_rate") = 0.0;
    parameters.doubles("apoptosis_resistance") = 1.0;
    parameters.doubles("ecm_production_rate_base") = 0.0;
    parameters.doubles("ecm_production_rate_boosted") = 0.0;
    parameters.doubles("mmp2_degradation_rate") = 0.0;
    parameters.doubles("mechanical_compaction_strength") = 0.0;
    parameters.doubles("compaction_ecm_increment") = 0.0;

    reset_all_fields(0.25, 0.0, 0.0, 0.0, 0.0, 0.6);

    const double cx = 0.0;
    const double cy = 0.0;
    place_tumor_cluster(pTumor, 5, cx, cy, 20.0);

    const double two_pi = 6.283185307179586;
    std::vector<Cell*> pscs;
    for (int i = 0; i < 10; ++i)
    {
        const double theta = two_pi * (static_cast<double>(i) / 10.0);
        const double radius = (i % 2 == 0) ? 30.0 : 40.0;
        Cell* p = create_cell(*pStroma);
        p->assign_position(std::vector<double>{cx + radius * std::cos(theta), cy + radius * std::sin(theta), 0.0});
        p->phenotype.motility.is_motile = false;
        p->custom_data[custom_index(p, "acta2_active")] = 0.0;
        p->custom_data[custom_index(p, "gli1_active")] = 0.0;
        pscs.push_back(p);
    }

    double t = 0.0;
    const double dt = 6.0;

    std::vector<double> previous_acta2(pscs.size(), 0.0);
    int activated_transitions = 0;
    int activated_with_gli1_at_transition = 0;

    for (int step = 0; step < 100; ++step)
    {
        advance_steps(1, dt, t);
        for (size_t i = 0; i < pscs.size(); ++i)
        {
            Cell* p = pscs[i];
            if (p == NULL || p->phenotype.death.dead) continue;
            const double act = p->custom_data[custom_index(p, "acta2_active")];
            const double gli = p->custom_data[custom_index(p, "gli1_active")];
            if (previous_acta2[i] == 0.0 && act == 1.0)
            {
                ++activated_transitions;
                if (gli > 0.5) ++activated_with_gli1_at_transition;
            }
            previous_acta2[i] = act;
        }
    }

    int activated = 0;
    for (size_t i = 0; i < pscs.size(); ++i)
    {
        Cell* p = pscs[i];
        if (p == NULL || p->phenotype.death.dead) continue;
        const double act = p->custom_data[custom_index(p, "acta2_active")];
        if (act == 1.0)
        {
            ++activated;
        }
    }

    const double shh_peak_before_off = max_field_value(shh_index);

    // Turn off SHH secretion, clear SHH field, then verify GLI1 reverses while ACTA2 stays ON.
    parameters.doubles("shh_secretion_rate") = 0.0;
    const int n_vox = static_cast<int>(microenvironment.number_of_voxels());
    for (int n = 0; n < n_vox; ++n)
    {
        microenvironment.density_vector(n)[shh_index] = 0.0;
    }

    advance_steps(20, dt, t);

    int activated_after = 0;
    int activated_gli1_off_after = 0;
    int caf_shh_nonzero_count = 0;

    for (size_t i = 0; i < pscs.size(); ++i)
    {
        Cell* p = pscs[i];
        if (p == NULL || p->phenotype.death.dead) continue;
        const double act = p->custom_data[custom_index(p, "acta2_active")];
        const double gli = p->custom_data[custom_index(p, "gli1_active")];
        if (act == 1.0)
        {
            ++activated_after;
            if (gli == 0.0) ++activated_gli1_off_after;

            const double caf_shh_rate = p->phenotype.secretion.secretion_rates[shh_index];
            if (caf_shh_rate > 0.0) ++caf_shh_nonzero_count;
        }
    }

    const double tgfb_peak = max_field_value(tgfb_index);
    const double shh_peak_after_off = max_field_value(shh_index);

    const bool pass_shh_sufficient = (activated > 0) && (tgfb_peak <= 1e-6);
    const bool pass_gli1_on_with_shh =
        (activated_transitions > 0) &&
        (activated_with_gli1_at_transition == activated_transitions);
    const bool pass_gli1_reversible = (activated_after > 0) && (activated_gli1_off_after == activated_after);
    const bool pass_acta2_irreversible = (activated_after == activated);
    const bool pass_caf_no_shh =
        (caf_shh_nonzero_count == 0) &&
        (shh_peak_after_off < shh_peak_before_off) &&
        (shh_peak_after_off <= 1e-3);

    std::cout << "E07 measurements"
              << " activated_step100=" << activated
              << " activated_transitions=" << activated_transitions
              << " activated_with_gli1_at_transition=" << activated_with_gli1_at_transition
              << " activated_after_shh_off=" << activated_after
              << " activated_gli1_off_after=" << activated_gli1_off_after
              << " tgfb_peak=" << tgfb_peak
              << " shh_peak_before_off=" << shh_peak_before_off
              << " shh_peak_after_off=" << shh_peak_after_off
              << " caf_nonzero_shh_rate_count=" << caf_shh_nonzero_count
              << std::endl;

    std::cout << "E07 checks"
              << " shh_alone_sufficient=" << (pass_shh_sufficient ? 1 : 0)
              << " gli1_on_when_shh_present=" << (pass_gli1_on_with_shh ? 1 : 0)
              << " gli1_reverts_when_shh_removed=" << (pass_gli1_reversible ? 1 : 0)
              << " acta2_stays_on=" << (pass_acta2_irreversible ? 1 : 0)
              << " cafs_never_secrete_shh=" << (pass_caf_no_shh ? 1 : 0)
              << std::endl;

    if (!(pass_shh_sufficient &&
          pass_gli1_on_with_shh &&
          pass_gli1_reversible &&
          pass_acta2_irreversible &&
          pass_caf_no_shh))
    {
        std::cout << "FAIL E07 integration" << std::endl;
        return 1;
    }

    std::cout << "PASS E07 integration" << std::endl;
    return 0;
}
