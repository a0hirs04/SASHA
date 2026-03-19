#include <cassert>
#include <iostream>

#include "feedback_loop_common.h"

using namespace feedback_loop_common;

namespace
{

struct Loop5Snapshot
{
    int step = 0;
    int live_tumor_n = 0;
    int live_exposed_n = 0;
    double abcb1_mean_all = 0.0;
    double abcb1_mean_exposed = 0.0;
    double abcb1_fraction_all = 0.0;
    double abcb1_fraction_exposed = 0.0;
    double intracellular_drug_mean_exposed = 0.0;
    double tgfb_secretion_mean_exposed = 0.0;
    int caf_count = 0;
};

Loop5Snapshot measure_snapshot(int step, int tumor_type, int stromal_type)
{
    Loop5Snapshot s;
    s.step = step;

    int abcb1_on_all = 0;
    int abcb1_on_exposed = 0;
    double abcb1_sum_all = 0.0;
    double abcb1_sum_exposed = 0.0;
    double drug_sum_exposed = 0.0;
    double tgfb_sec_sum_exposed = 0.0;
    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->phenotype.death.dead) continue;
        if (pCell->type != tumor_type) continue;

        const int abcb1_idx = custom_index(pCell, "abcb1_active");
        const int intra_idx = custom_index(pCell, "intracellular_drug");
        const int texp_idx = custom_index(pCell, "time_since_drug_exposure");
        const bool exposed = (texp_idx >= 0 && pCell->custom_data[texp_idx] >= 0.0);

        if (abcb1_idx >= 0)
        {
            const double abcb1_value = pCell->custom_data[abcb1_idx];
            abcb1_sum_all += abcb1_value;
            if (abcb1_value > 0.3) ++abcb1_on_all;
        }
        if (exposed)
        {
            ++s.live_exposed_n;
            if (abcb1_idx >= 0)
            {
                const double abcb1_value = pCell->custom_data[abcb1_idx];
                abcb1_sum_exposed += abcb1_value;
                if (abcb1_value > 0.3) ++abcb1_on_exposed;
            }
            if (intra_idx >= 0) drug_sum_exposed += pCell->custom_data[intra_idx];
            if (tgfb_index >= 0 &&
                tgfb_index < static_cast<int>(pCell->phenotype.secretion.secretion_rates.size()))
            {
                tgfb_sec_sum_exposed += pCell->phenotype.secretion.secretion_rates[tgfb_index];
            }
        }
        ++s.live_tumor_n;
    }

    if (s.live_tumor_n > 0)
    {
        s.abcb1_mean_all = abcb1_sum_all / static_cast<double>(s.live_tumor_n);
        s.abcb1_fraction_all = static_cast<double>(abcb1_on_all) / static_cast<double>(s.live_tumor_n);
    }
    if (s.live_exposed_n > 0)
    {
        s.abcb1_mean_exposed =
            abcb1_sum_exposed / static_cast<double>(s.live_exposed_n);
        s.abcb1_fraction_exposed =
            static_cast<double>(abcb1_on_exposed) / static_cast<double>(s.live_exposed_n);
        s.intracellular_drug_mean_exposed =
            drug_sum_exposed / static_cast<double>(s.live_exposed_n);
        s.tgfb_secretion_mean_exposed =
            tgfb_sec_sum_exposed / static_cast<double>(s.live_exposed_n);
    }

    s.caf_count = count_activated_cafs(stromal_type);
    return s;
}

} // namespace

int main()
{
    initialize_world();

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    Cell_Definition* pStroma = find_cell_definition("stromal_cell");
    assert(pTumor != NULL);
    assert(pStroma != NULL);
    const int tumor_type = pTumor->type;
    const int stromal_type = pStroma->type;

    // Build and maintain stromal barrier.
    parameters.doubles("tgfb_secretion_rate") = 0.45;
    parameters.doubles("shh_secretion_rate") = 0.35;
    parameters.doubles("hif1a_tgfb_amplification_factor") = 1.4;
    parameters.doubles("caf_tgfb_secretion_rate") = 0.30;
    parameters.doubles("tgfb_activation_threshold") = 0.18;
    parameters.doubles("shh_activation_threshold") = 0.05;
    parameters.doubles("ecm_production_rate_base") = 0.004;
    parameters.doubles("ecm_production_rate_boosted") = 0.007;
    parameters.doubles("ecm_ha_fraction_default") = 0.6;
    parameters.doubles("mmp2_degradation_rate") = 0.0;
    parameters.doubles("emt_induction_threshold") = 10.0;

    // Keep population stable so adaptation dominates readouts.
    parameters.doubles("base_proliferation_rate") = 0.0;
    parameters.doubles("caf_proliferation_rate") = 0.0;
    parameters.doubles("psc_proliferation_rate") = 0.0;
    parameters.doubles("apoptosis_resistance") = 0.98;

    // Drug-response module settings.
    parameters.doubles("drug_uptake_rate") = 0.006;
    parameters.doubles("drug_stress_threshold") = 0.2;
    parameters.doubles("efflux_induction_delay") = 80.0;
    parameters.doubles("efflux_strength") = 0.05;
    parameters.doubles("nrf2_decay_rate") = 0.01;
    parameters.doubles("drug_natural_decay_rate") = 0.0;
    parameters.doubles("hif1a_nrf2_priming_bonus") = 0.0;

    // Drug kill in Module 4.
    parameters.doubles("drug_kill_coefficient") = 0.0015;
    parameters.doubles("efflux_drug_reduction") = 0.7;
    microenvironment.diffusion_coefficients[drug_index] = 500.0;

    reset_all_fields(0.2, 0.0, 0.0, 0.0, 0.0, 0.6);

    const double cx = 0.0;
    const double cy = 0.0;
    // Treated cohort near center + untreated distant cohort for delayed population-level adaptation.
    place_tumor_cluster(pTumor, 20, cx, cy, 60.0);
    place_tumor_cluster(pTumor, 20, -500.0, 0.0, 40.0);
    place_stroma_ring(pStroma, 40, cx, cy, 130.0, 220.0);

    double t = 0.0;
    const double dt = 1.0;

    // Mature barrier before treatment.
    advance_steps(150, dt, t);

    Loop5Snapshot s50;
    Loop5Snapshot s150;
    Loop5Snapshot s300;

    // Treatment window: apply boundary drug continuously.
    for (int step = 1; step <= 300; ++step)
    {
        set_annulus_field(drug_index, cx, cy, 180.0, 300.0, 1.0);
        advance_steps(1, dt, t);

        if (step == 50) s50 = measure_snapshot(step, tumor_type, stromal_type);
        if (step == 150) s150 = measure_snapshot(step, tumor_type, stromal_type);
        if (step == 300) s300 = measure_snapshot(step, tumor_type, stromal_type);
    }

    const bool pass =
        (s50.live_tumor_n > 0 && s150.live_tumor_n > 0 && s300.live_tumor_n > 0) &&
        (s50.live_exposed_n > 0 && s300.live_exposed_n > 0) &&
        (s50.abcb1_fraction_all < 0.2) &&
        (s50.abcb1_mean_all < s150.abcb1_mean_all) &&
        (s150.abcb1_mean_all < s300.abcb1_mean_all) &&
        (s300.abcb1_fraction_exposed > 0.3) &&
        (s300.intracellular_drug_mean_exposed < s150.intracellular_drug_mean_exposed) &&
        (s300.tgfb_secretion_mean_exposed > 0.0) &&
        (s300.caf_count >= s50.caf_count);

    std::cout << "LOOP5 step50 live_tumor=" << s50.live_tumor_n
              << " exposed_n=" << s50.live_exposed_n
              << " abcb1_mean_all=" << s50.abcb1_mean_all
              << " abcb1_mean_exposed=" << s50.abcb1_mean_exposed
              << " abcb1_fraction_all=" << s50.abcb1_fraction_all
              << " abcb1_fraction_exposed=" << s50.abcb1_fraction_exposed
              << " intracellular_drug_mean_exposed=" << s50.intracellular_drug_mean_exposed
              << " tgfb_secretion_mean_exposed=" << s50.tgfb_secretion_mean_exposed
              << " caf_count=" << s50.caf_count << std::endl;
    std::cout << "LOOP5 step150 live_tumor=" << s150.live_tumor_n
              << " exposed_n=" << s150.live_exposed_n
              << " abcb1_mean_all=" << s150.abcb1_mean_all
              << " abcb1_mean_exposed=" << s150.abcb1_mean_exposed
              << " abcb1_fraction_all=" << s150.abcb1_fraction_all
              << " abcb1_fraction_exposed=" << s150.abcb1_fraction_exposed
              << " intracellular_drug_mean_exposed=" << s150.intracellular_drug_mean_exposed
              << " tgfb_secretion_mean_exposed=" << s150.tgfb_secretion_mean_exposed
              << " caf_count=" << s150.caf_count << std::endl;
    std::cout << "LOOP5 step300 live_tumor=" << s300.live_tumor_n
              << " exposed_n=" << s300.live_exposed_n
              << " abcb1_mean_all=" << s300.abcb1_mean_all
              << " abcb1_mean_exposed=" << s300.abcb1_mean_exposed
              << " abcb1_fraction_all=" << s300.abcb1_fraction_all
              << " abcb1_fraction_exposed=" << s300.abcb1_fraction_exposed
              << " intracellular_drug_mean_exposed=" << s300.intracellular_drug_mean_exposed
              << " tgfb_secretion_mean_exposed=" << s300.tgfb_secretion_mean_exposed
              << " caf_count=" << s300.caf_count << std::endl;

    std::cout << "LOOP5 " << (pass ? "PASS" : "FAIL") << std::endl;
    if (!pass) return 1;

    std::cout << "PASS loop5_drug_resistance_test" << std::endl;
    return 0;
}
