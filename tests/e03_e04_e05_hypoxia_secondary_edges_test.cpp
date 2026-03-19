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

bool nearly_equal(double a, double b, double eps = 1e-12)
{
    return std::fabs(a - b) <= eps;
}

} // namespace

int main()
{
    initialize_world();

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    assert(pTumor != NULL);

    // E03 + E04 setup: compare HIF1A OFF vs ON with same drug and all else equal.
    parameters.doubles("base_proliferation_rate") = 0.01;
    parameters.doubles("go_grow_penalty") = 0.0;
    parameters.doubles("hypoxia_proliferation_modifier") = 0.3;
    parameters.doubles("contact_inhibition_threshold") = 1e9;
    parameters.doubles("apoptosis_resistance") = 0.85;
    parameters.doubles("hypoxia_death_resistance_bonus") = 0.5;
    parameters.doubles("drug_kill_coefficient") = 0.0001;
    parameters.doubles("efflux_drug_reduction") = 0.0;

    // E05 setup.
    parameters.doubles("drug_uptake_rate") = 0.0;
    parameters.doubles("drug_stress_threshold") = 0.03;
    parameters.doubles("efflux_induction_delay") = 1e9;
    parameters.doubles("efflux_strength") = 0.0;
    parameters.doubles("nrf2_decay_rate") = 0.01;
    parameters.doubles("drug_natural_decay_rate") = 0.0;
    parameters.doubles("hif1a_nrf2_priming_bonus") = 0.02;

    parameters.doubles("tgfb_secretion_rate") = 0.0;
    parameters.doubles("shh_secretion_rate") = 0.0;
    parameters.doubles("ecm_production_rate_base") = 0.0;
    parameters.doubles("ecm_production_rate_boosted") = 0.0;
    parameters.doubles("mmp2_degradation_rate") = 0.0;
    parameters.doubles("mechanical_compaction_strength") = 0.0;
    parameters.doubles("compaction_ecm_increment") = 0.0;
    parameters.doubles("emt_induction_threshold") = 10.0;

    reset_all_fields(0.25, 0.0, 0.0, 0.0, 0.0, 0.6);

    // E03/E04 pair.
    Cell* hif_off = create_cell(*pTumor);
    hif_off->assign_position(std::vector<double>{-100.0, 0.0, 0.0});
    Cell* hif_on = create_cell(*pTumor);
    hif_on->assign_position(std::vector<double>{100.0, 0.0, 0.0});

    const int hif_idx_off = custom_index(hif_off, "hif1a_active");
    const int hif_idx_on = custom_index(hif_on, "hif1a_active");
    const int zeb_off = custom_index(hif_off, "zeb1_active");
    const int zeb_on = custom_index(hif_on, "zeb1_active");
    const int abcb1_off = custom_index(hif_off, "abcb1_active");
    const int abcb1_on = custom_index(hif_on, "abcb1_active");
    const int intra_off = custom_index(hif_off, "intracellular_drug");
    const int intra_on = custom_index(hif_on, "intracellular_drug");
    const int pressure_off = custom_index(hif_off, "mechanical_pressure");
    const int pressure_on = custom_index(hif_on, "mechanical_pressure");

    hif_off->custom_data[hif_idx_off] = 0.0;
    hif_on->custom_data[hif_idx_on] = 1.0;
    hif_off->custom_data[zeb_off] = 0.0;
    hif_on->custom_data[zeb_on] = 0.0;
    hif_off->custom_data[abcb1_off] = 0.0;
    hif_on->custom_data[abcb1_on] = 0.0;
    hif_off->custom_data[intra_off] = 0.1;
    hif_on->custom_data[intra_on] = 0.1;
    hif_off->custom_data[pressure_off] = 0.0;
    hif_on->custom_data[pressure_on] = 0.0;

    module4_proliferation_death(hif_off, hif_off->phenotype, 1.0, ModulePhase::DECISION);
    module4_proliferation_death(hif_on, hif_on->phenotype, 1.0, ModulePhase::DECISION);

    const int apop_idx_off = hif_off->phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);
    const int apop_idx_on = hif_on->phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);

    const double prolif_off = hif_off->phenotype.cycle.data.transition_rate(0, 0);
    const double prolif_on = hif_on->phenotype.cycle.data.transition_rate(0, 0);
    const double apop_off = hif_off->phenotype.death.rates[apop_idx_off];
    const double apop_on = hif_on->phenotype.death.rates[apop_idx_on];

    const bool pass_e03 = (prolif_on < prolif_off);
    const bool pass_e04 = (apop_on < apop_off);
    const bool pass_e03e04_simultaneous = pass_e03 && pass_e04;

    // E05: NRF2 priming by HIF1A near threshold.
    Cell* e05_a = create_cell(*pTumor);
    e05_a->assign_position(std::vector<double>{-200.0, 0.0, 0.0});
    Cell* e05_b = create_cell(*pTumor);
    e05_b->assign_position(std::vector<double>{200.0, 0.0, 0.0});

    const int hif_a = custom_index(e05_a, "hif1a_active");
    const int hif_b = custom_index(e05_b, "hif1a_active");
    const int nrf2_a = custom_index(e05_a, "nrf2_active");
    const int nrf2_b = custom_index(e05_b, "nrf2_active");
    const int intra_a = custom_index(e05_a, "intracellular_drug");
    const int intra_b = custom_index(e05_b, "intracellular_drug");
    const int texp_a = custom_index(e05_a, "time_since_drug_exposure");
    const int texp_b = custom_index(e05_b, "time_since_drug_exposure");

    e05_a->custom_data[hif_a] = 0.0;
    e05_b->custom_data[hif_b] = 1.0;
    e05_a->custom_data[nrf2_a] = 0.0;
    e05_b->custom_data[nrf2_b] = 0.0;
    e05_a->custom_data[intra_a] = 0.02;
    e05_b->custom_data[intra_b] = 0.02;
    e05_a->custom_data[texp_a] = -1.0;
    e05_b->custom_data[texp_b] = -1.0;

    e05_a->nearest_density_vector()[drug_index] = 0.0;
    e05_b->nearest_density_vector()[drug_index] = 0.0;

    module7_drug_response(e05_a, e05_a->phenotype, 1.0, ModulePhase::SENSING);
    module7_drug_response(e05_b, e05_b->phenotype, 1.0, ModulePhase::SENSING);

    const bool pass_e05_positive =
        nearly_equal(e05_a->custom_data[nrf2_a], 0.0) &&
        (e05_b->custom_data[nrf2_b] > e05_a->custom_data[nrf2_a]);

    // Mandatory negative: with no intracellular drug, HIF1A alone must not activate NRF2.
    Cell* e05_neg = create_cell(*pTumor);
    e05_neg->assign_position(std::vector<double>{300.0, 0.0, 0.0});
    const int hif_neg = custom_index(e05_neg, "hif1a_active");
    const int nrf2_neg = custom_index(e05_neg, "nrf2_active");
    const int intra_neg = custom_index(e05_neg, "intracellular_drug");
    const int texp_neg = custom_index(e05_neg, "time_since_drug_exposure");

    e05_neg->custom_data[hif_neg] = 1.0;
    e05_neg->custom_data[nrf2_neg] = 0.0;
    e05_neg->custom_data[intra_neg] = 0.0;
    e05_neg->custom_data[texp_neg] = -1.0;
    e05_neg->nearest_density_vector()[drug_index] = 0.0;

    module7_drug_response(e05_neg, e05_neg->phenotype, 1.0, ModulePhase::SENSING);
    const bool pass_e05_negative = nearly_equal(e05_neg->custom_data[nrf2_neg], 0.0);

    const bool pass =
        pass_e03 &&
        pass_e04 &&
        pass_e03e04_simultaneous &&
        pass_e05_positive &&
        pass_e05_negative;

    std::cout << "P5A measurements"
              << " prolif_hif_off=" << prolif_off
              << " prolif_hif_on=" << prolif_on
              << " apop_hif_off=" << apop_off
              << " apop_hif_on=" << apop_on
              << " nrf2_e05_a=" << e05_a->custom_data[nrf2_a]
              << " nrf2_e05_b=" << e05_b->custom_data[nrf2_b]
              << " nrf2_e05_neg=" << e05_neg->custom_data[nrf2_neg]
              << std::endl;

    std::cout << "P5A checks"
              << " E03_prolif_hif_on_lt_off=" << (pass_e03 ? 1 : 0)
              << " E04_apop_hif_on_lt_off=" << (pass_e04 ? 1 : 0)
              << " E03E04_simultaneous=" << (pass_e03e04_simultaneous ? 1 : 0)
              << " E05_priming_positive=" << (pass_e05_positive ? 1 : 0)
              << " E05_hif_only_negative=" << (pass_e05_negative ? 1 : 0)
              << std::endl;

    if (!pass)
    {
        std::cout << "FAIL P5 Cluster A" << std::endl;
        return 1;
    }

    std::cout << "PASS P5 Cluster A" << std::endl;
    return 0;
}
