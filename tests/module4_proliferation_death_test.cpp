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

int get_apoptosis_index(Phenotype& phenotype)
{
    return phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);
}

} // namespace

int main()
{
    const bool xml_status = load_PhysiCell_config_file("config/PhysiCell_settings.xml");
    assert(xml_status);

    setup_microenvironment();
    create_cell_container_for_microenvironment(microenvironment, 30.0);
    create_cell_types();

    parameters.doubles("base_proliferation_rate") = 0.01;
    parameters.doubles("apoptosis_resistance") = 0.85;
    parameters.doubles("go_grow_penalty") = 0.3;
    parameters.doubles("hypoxia_proliferation_modifier") = 0.2;
    parameters.doubles("hypoxia_death_resistance_bonus") = 0.1;
    parameters.doubles("contact_inhibition_threshold") = 5.0;
    parameters.doubles("drug_kill_coefficient") = 0.01;
    parameters.doubles("efflux_drug_reduction") = 0.7;
    parameters.doubles("caf_proliferation_rate") = 0.003;
    parameters.doubles("gli1_proliferation_boost") = 1.3;
    parameters.doubles("psc_proliferation_rate") = 0.0;

    const double baseline_death = (1.0 - 0.85) * 0.001; // 0.00015

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    Cell_Definition* pStroma = find_cell_definition("stromal_cell");
    assert(pTumor != NULL);
    assert(pStroma != NULL);

    Cell* tumor = create_cell(*pTumor);
    tumor->assign_position(std::vector<double>{100.0, 100.0, 0.0});

    const int tumor_apop_idx = get_apoptosis_index(tumor->phenotype);
    assert(tumor_apop_idx >= 0);
    const int t_hif1a = tumor->custom_data.find_variable_index("hif1a_active");
    const int t_zeb1 = tumor->custom_data.find_variable_index("zeb1_active");
    const int t_abcb1 = tumor->custom_data.find_variable_index("abcb1_active");
    const int t_drug = tumor->custom_data.find_variable_index("intracellular_drug");
    const int t_pressure = tumor->custom_data.find_variable_index("mechanical_pressure");
    assert(t_hif1a >= 0);
    assert(t_zeb1 >= 0);
    assert(t_abcb1 >= 0);
    assert(t_drug >= 0);
    assert(t_pressure >= 0);

    // Test A — Tumor baseline.
    tumor->custom_data[t_hif1a] = 0.0;
    tumor->custom_data[t_zeb1] = 0.0;
    tumor->custom_data[t_abcb1] = 0.0;
    tumor->custom_data[t_drug] = 0.0;
    tumor->custom_data[t_pressure] = 0.0;
    module4_proliferation_death(tumor, tumor->phenotype, 1.0, ModulePhase::DECISION);
    assert(nearly_equal(tumor->phenotype.cycle.data.transition_rate(0, 0), 0.01));
    assert(nearly_equal(tumor->phenotype.death.rates[tumor_apop_idx], baseline_death));
    std::cout << "PASS Test A" << std::endl;

    // Test B — Tumor mechanical compression.
    tumor->custom_data[t_pressure] = 6.0;
    module4_proliferation_death(tumor, tumor->phenotype, 1.0, ModulePhase::DECISION);
    assert(nearly_equal(tumor->phenotype.cycle.data.transition_rate(0, 0), 0.0));
    assert(nearly_equal(tumor->phenotype.death.rates[tumor_apop_idx], baseline_death));
    std::cout << "PASS Test B" << std::endl;

    // Test C — Tumor EMT penalty.
    tumor->custom_data[t_pressure] = 0.0;
    tumor->custom_data[t_zeb1] = 1.0;
    tumor->custom_data[t_hif1a] = 0.0;
    tumor->custom_data[t_drug] = 0.0;
    tumor->custom_data[t_abcb1] = 0.0;
    module4_proliferation_death(tumor, tumor->phenotype, 1.0, ModulePhase::DECISION);
    assert(nearly_equal(tumor->phenotype.cycle.data.transition_rate(0, 0), 0.007));
    std::cout << "PASS Test C" << std::endl;

    // Test D — Tumor drug kill, no efflux.
    tumor->custom_data[t_zeb1] = 0.0;
    tumor->custom_data[t_hif1a] = 0.0;
    tumor->custom_data[t_abcb1] = 0.0;
    tumor->custom_data[t_drug] = 0.5;
    module4_proliferation_death(tumor, tumor->phenotype, 1.0, ModulePhase::DECISION);
    assert(nearly_equal(tumor->phenotype.death.rates[tumor_apop_idx], baseline_death + 0.005));
    std::cout << "PASS Test D" << std::endl;

    // Test E — Tumor drug kill with ABCB1 efflux.
    tumor->custom_data[t_abcb1] = 1.0;
    tumor->custom_data[t_drug] = 0.5;
    module4_proliferation_death(tumor, tumor->phenotype, 1.0, ModulePhase::DECISION);
    assert(nearly_equal(tumor->phenotype.death.rates[tumor_apop_idx], baseline_death + 0.0015));
    std::cout << "PASS Test E" << std::endl;

    // Test F — Activated CAF with GLI1.
    Cell* stroma_f = create_cell(*pStroma);
    stroma_f->assign_position(std::vector<double>{140.0, 100.0, 0.0});
    const int s_apop_idx_f = get_apoptosis_index(stroma_f->phenotype);
    const int s_acta2_f = stroma_f->custom_data.find_variable_index("acta2_active");
    const int s_gli1_f = stroma_f->custom_data.find_variable_index("gli1_active");
    assert(s_apop_idx_f >= 0);
    assert(s_acta2_f >= 0);
    assert(s_gli1_f >= 0);
    stroma_f->custom_data[s_acta2_f] = 1.0;
    stroma_f->custom_data[s_gli1_f] = 1.0;
    module4_proliferation_death(stroma_f, stroma_f->phenotype, 1.0, ModulePhase::DECISION);
    assert(nearly_equal(stroma_f->phenotype.cycle.data.transition_rate(0, 0), 0.0039));
    assert(nearly_equal(stroma_f->phenotype.death.rates[s_apop_idx_f], 0.0));
    std::cout << "PASS Test F" << std::endl;

    // Test G — Stromal cell ignores drug.
    Cell* stroma_g = create_cell(*pStroma);
    stroma_g->assign_position(std::vector<double>{180.0, 100.0, 0.0});
    const int s_apop_idx_g = get_apoptosis_index(stroma_g->phenotype);
    const int s_acta2_g = stroma_g->custom_data.find_variable_index("acta2_active");
    const int s_drug_g = stroma_g->custom_data.find_variable_index("intracellular_drug");
    assert(s_apop_idx_g >= 0);
    assert(s_acta2_g >= 0);
    stroma_g->custom_data[s_acta2_g] = 1.0;
    if (s_drug_g >= 0)
    {
        stroma_g->custom_data[s_drug_g] = 0.9;
    }
    module4_proliferation_death(stroma_g, stroma_g->phenotype, 1.0, ModulePhase::DECISION);
    assert(nearly_equal(stroma_g->phenotype.death.rates[s_apop_idx_g], 0.0));
    std::cout << "PASS Test G" << std::endl;

    std::cout << "PASS module4_proliferation_death_test" << std::endl;
    return 0;
}
