#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
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

void set_local_density(Cell* pCell, int density_index, double value)
{
    const int voxel = voxel_index_for_cell(pCell);
    microenvironment.density_vector(voxel)[density_index] = value;
}

double get_density_at_voxel(int voxel, int density_index)
{
    return microenvironment.density_vector(voxel)[density_index];
}

std::string phase_to_string(ModulePhase phase)
{
    switch (phase)
    {
        case ModulePhase::SENSING: return "SENSING";
        case ModulePhase::DECISION: return "DECISION";
        case ModulePhase::WRITE: return "WRITE";
        default: return "UNKNOWN";
    }
}

std::string join(const std::vector<std::string>& values)
{
    if (values.empty()) return "none";
    std::string out;
    for (size_t i = 0; i < values.size(); ++i)
    {
        if (i > 0) out += ",";
        out += values[i];
    }
    return out;
}

bool same_strings(std::vector<std::string> a, std::vector<std::string> b)
{
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    return a == b;
}

} // namespace

int main()
{
    const bool xml_status = load_PhysiCell_config_file("config/PhysiCell_settings.xml");
    assert(xml_status);

    setup_microenvironment();
    create_cell_container_for_microenvironment(microenvironment, 30.0);
    create_cell_types();

    assert(oxygen_index >= 0);
    assert(tgfb_index >= 0);
    assert(shh_index >= 0);
    assert(ecm_index >= 0);
    assert(drug_index >= 0);

    parameters.doubles("hypoxia_response_threshold") = 0.35;

    parameters.doubles("tgfb_activation_threshold") = 0.6;
    parameters.doubles("shh_activation_threshold") = 0.1;

    parameters.doubles("drug_uptake_rate") = 0.0;
    parameters.doubles("drug_stress_threshold") = 0.03;
    parameters.doubles("efflux_induction_delay") = 0.0;
    parameters.doubles("efflux_strength") = 0.0;
    parameters.doubles("nrf2_decay_rate") = 0.01;
    parameters.doubles("drug_natural_decay_rate") = 0.0;
    parameters.doubles("hif1a_nrf2_priming_bonus") = 0.0;

    parameters.doubles("emt_induction_threshold") = 0.4;
    parameters.doubles("emt_activation_delay") = 0.0;
    parameters.ints("emt_phenotype_extent") = 2;
    parameters.doubles("hif1a_emt_boost") = 0.2;
    parameters.doubles("motility_epithelial") = 0.1;
    parameters.doubles("motility_mesenchymal_low") = 0.5;
    parameters.doubles("motility_mesenchymal_med") = 1.0;
    parameters.doubles("motility_mesenchymal_high") = 2.0;
    parameters.doubles("adhesion_epithelial") = 5.0;
    parameters.doubles("adhesion_mesenchymal_low") = 3.0;
    parameters.doubles("adhesion_mesenchymal_med") = 1.5;
    parameters.doubles("adhesion_mesenchymal_high") = 0.5;

    parameters.doubles("base_proliferation_rate") = 0.01;
    parameters.doubles("apoptosis_resistance") = 0.85;
    parameters.doubles("go_grow_penalty") = 0.3;
    parameters.doubles("hypoxia_proliferation_modifier") = 0.0;
    parameters.doubles("hypoxia_death_resistance_bonus") = 0.0;
    parameters.doubles("contact_inhibition_threshold") = 999.0;
    parameters.doubles("drug_kill_coefficient") = 0.01;
    parameters.doubles("efflux_drug_reduction") = 0.7;
    parameters.doubles("caf_proliferation_rate") = 0.003;
    parameters.doubles("gli1_proliferation_boost") = 1.3;
    parameters.doubles("psc_proliferation_rate") = 0.0;

    parameters.doubles("tgfb_secretion_rate") = 0.5;
    parameters.doubles("shh_secretion_rate") = 0.3;
    parameters.doubles("hif1a_tgfb_amplification_factor") = 1.5;
    parameters.doubles("caf_tgfb_secretion_rate") = 0.3;

    parameters.doubles("ecm_production_rate_base") = 0.02;
    parameters.doubles("ecm_production_rate_boosted") = 0.03;
    parameters.doubles("ecm_ha_fraction_default") = 0.6;
    parameters.doubles("mmp2_degradation_rate") = 0.01;
    parameters.doubles("ecm_natural_decay_rate") = 0.0;

    parameters.doubles("mechanical_compaction_strength") = 0.0;
    parameters.doubles("compaction_ecm_increment") = 0.0;

    const double dt = 1.0;

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    Cell_Definition* pStroma = find_cell_definition("stromal_cell");
    assert(pTumor != NULL);
    assert(pStroma != NULL);

    Cell* tumor = create_cell(*pTumor);
    tumor->assign_position(std::vector<double>{100.0, 100.0, 0.0});

    Cell* psc = create_cell(*pStroma);
    psc->assign_position(std::vector<double>{220.0, 100.0, 0.0});

    const int tumor_hif1a = tumor->custom_data.find_variable_index("hif1a_active");
    const int tumor_zeb1 = tumor->custom_data.find_variable_index("zeb1_active");
    const int tumor_mmp2 = tumor->custom_data.find_variable_index("mmp2_active");
    const int tumor_nrf2 = tumor->custom_data.find_variable_index("nrf2_active");
    const int tumor_abcb1 = tumor->custom_data.find_variable_index("abcb1_active");
    const int tumor_intra = tumor->custom_data.find_variable_index("intracellular_drug");
    const int tumor_texp = tumor->custom_data.find_variable_index("time_since_drug_exposure");
    assert(tumor_hif1a >= 0);
    assert(tumor_zeb1 >= 0);
    assert(tumor_mmp2 >= 0);
    assert(tumor_nrf2 >= 0);
    assert(tumor_abcb1 >= 0);
    assert(tumor_intra >= 0);
    assert(tumor_texp >= 0);

    const int psc_acta2 = psc->custom_data.find_variable_index("acta2_active");
    const int psc_gli1 = psc->custom_data.find_variable_index("gli1_active");
    const int psc_ecm_rate = psc->custom_data.find_variable_index("ecm_production_rate");
    assert(psc_acta2 >= 0);
    assert(psc_gli1 >= 0);
    assert(psc_ecm_rate >= 0);

    const int tumor_voxel = voxel_index_for_cell(tumor);
    const int psc_voxel = voxel_index_for_cell(psc);

    set_local_density(tumor, oxygen_index, 0.02);
    set_local_density(tumor, tgfb_index, 0.25);
    set_local_density(tumor, shh_index, 0.0);
    set_local_density(tumor, drug_index, 0.5);
    microenvironment.density_vector(tumor_voxel)[ecm_index] = 0.5;
    set_ecm_ha_fraction(tumor_voxel, 0.4);

    set_local_density(psc, oxygen_index, 0.08);
    set_local_density(psc, tgfb_index, 0.5);
    set_local_density(psc, shh_index, 0.3);
    set_local_density(psc, drug_index, 0.0);
    microenvironment.density_vector(psc_voxel)[ecm_index] = 0.0;
    set_ecm_ha_fraction(psc_voxel, 0.5);

    tumor->custom_data[tumor_hif1a] = 0.0;
    tumor->custom_data[tumor_zeb1] = 0.0;
    tumor->custom_data[tumor_mmp2] = 0.0;
    tumor->custom_data[tumor_nrf2] = 1.0;
    tumor->custom_data[tumor_abcb1] = 0.0;
    tumor->custom_data[tumor_intra] = 0.5;
    tumor->custom_data[tumor_texp] = 0.0;
    tumor->state.simple_pressure = 0.0;

    psc->custom_data[psc_acta2] = 0.0;
    psc->state.simple_pressure = 0.0;

    // Test 1 — Order dependency check.
    custom_function(tumor, tumor->phenotype, dt);
    custom_function(psc, psc->phenotype, dt);

    // (a) Module1 -> Module5 dependency.
    assert(nearly_equal(tumor->custom_data[tumor_hif1a], 1.0));
    assert(nearly_equal(tumor->custom_data[tumor_zeb1], 1.0));

    // (d) Module5 -> Module4 dependency (go-grow penalty applies immediately).
    assert(nearly_equal(tumor->phenotype.cycle.data.transition_rate(0, 0), 0.007));

    // (c) Module7 -> Module4 dependency (ABCB1 reduces effective drug in same step).
    assert(tumor->custom_data[tumor_abcb1] > 0.0);
    const int tumor_apop_index =
        tumor->phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);
    assert(tumor_apop_index >= 0);
    const double emt_death_increase = parameters.doubles("emt_death_increase");
    const double expected_apoptosis =
        (1.0 - 0.85) * 0.001 +
        (0.5 * (1.0 - 0.7 * tumor->custom_data[tumor_abcb1])) * 0.01 +
        emt_death_increase;
    assert(nearly_equal(tumor->phenotype.death.rates[tumor_apop_index], expected_apoptosis));

    // (e) Module5 -> Module6 dependency (MMP2 flag drives same-step ECM degradation).
    assert(nearly_equal(tumor->custom_data[tumor_mmp2], 1.0));
    assert(nearly_equal(get_density_at_voxel(tumor_voxel, ecm_index), 0.49));

    // (b) Module3 -> Module6 dependency (ACTA2 activation drives same-step ECM production).
    assert(nearly_equal(psc->custom_data[psc_acta2], 1.0));
    assert(psc->custom_data[psc_gli1] > 0.0);
    assert(psc->custom_data[psc_ecm_rate] > parameters.doubles("ecm_production_rate_base"));
    assert(nearly_equal(
        get_density_at_voxel(psc_voxel, ecm_index),
        psc->custom_data[psc_ecm_rate] * dt));
    std::cout << "PASS Test 1" << std::endl;

    // Test 2 — Phase violation and order trace check.
    Cell* tumor_trace = create_cell(*pTumor);
    tumor_trace->assign_position(std::vector<double>{340.0, 100.0, 0.0});
    const int trace_voxel = voxel_index_for_cell(tumor_trace);
    set_local_density(tumor_trace, oxygen_index, 0.02);
    set_local_density(tumor_trace, tgfb_index, 0.25);
    set_local_density(tumor_trace, shh_index, 0.0);
    set_local_density(tumor_trace, drug_index, 0.5);
    microenvironment.density_vector(trace_voxel)[ecm_index] = 0.5;
    set_ecm_ha_fraction(trace_voxel, 0.4);
    tumor_trace->custom_data[tumor_nrf2] = 1.0;
    tumor_trace->custom_data[tumor_abcb1] = 0.0;
    tumor_trace->custom_data[tumor_intra] = 0.5;
    tumor_trace->custom_data[tumor_texp] = 0.0;

    set_module_trace_enabled(true);
    clear_module_trace_log();
    custom_function(tumor_trace, tumor_trace->phenotype, dt);
    set_module_trace_enabled(false);

    const std::vector<ModuleTraceEntry> trace = get_module_trace_log();
    assert(trace.size() == 8);

    const std::vector<std::string> expected_modules = {
        "module1_oxygen_sensing",
        "module3_stromal_activation",
        "module7_drug_response",
        "module5_emt_engine",
        "module4_proliferation_death",
        "module2_paracrine_secretion",
        "module6_ecm_production",
        "module8_mechanical_compaction"
    };
    const std::vector<ModulePhase> expected_phases = {
        ModulePhase::SENSING,
        ModulePhase::SENSING,
        ModulePhase::SENSING,
        ModulePhase::DECISION,
        ModulePhase::DECISION,
        ModulePhase::WRITE,
        ModulePhase::WRITE,
        ModulePhase::WRITE
    };

    for (size_t i = 0; i < trace.size(); ++i)
    {
        assert(trace[i].module_name == expected_modules[i]);
        assert(trace[i].phase == expected_phases[i]);
    }

    for (size_t i = 0; i < trace.size(); ++i)
    {
        if (trace[i].phase == ModulePhase::SENSING || trace[i].phase == ModulePhase::DECISION)
        {
            assert(trace[i].field_writes.empty());
        }
    }

    for (size_t i = 0; i < trace.size(); ++i)
    {
        if (trace[i].module_name == "module2_paracrine_secretion")
        {
            assert(same_strings(trace[i].field_writes, std::vector<std::string>{"tgfb", "shh"}));
        }
        if (trace[i].module_name == "module6_ecm_production")
        {
            assert(same_strings(trace[i].field_writes, std::vector<std::string>{"ecm_density", "ecm_ha_fraction"}));
        }
        if (trace[i].module_name == "module8_mechanical_compaction")
        {
            assert(same_strings(trace[i].field_writes, std::vector<std::string>{"ecm_density"}));
        }
    }

    size_t first_decision = trace.size();
    size_t first_write = trace.size();
    size_t last_sensing = 0;
    size_t last_decision = 0;
    for (size_t i = 0; i < trace.size(); ++i)
    {
        if (trace[i].phase == ModulePhase::SENSING) last_sensing = i;
        if (trace[i].phase == ModulePhase::DECISION)
        {
            last_decision = i;
            if (first_decision == trace.size()) first_decision = i;
        }
        if (trace[i].phase == ModulePhase::WRITE && first_write == trace.size())
        {
            first_write = i;
        }
    }
    assert(last_sensing < first_decision);
    assert(last_decision < first_write);

    for (size_t i = 0; i < trace.size(); ++i)
    {
        std::cout << "LOG module=" << trace[i].module_name
                  << " phase=" << phase_to_string(trace[i].phase)
                  << " reads=" << join(trace[i].field_reads)
                  << " writes=" << join(trace[i].field_writes)
                  << std::endl;
    }

    std::cout << "PASS Test 2" << std::endl;
    std::cout << "PASS module_execution_order_test" << std::endl;
    return 0;
}
