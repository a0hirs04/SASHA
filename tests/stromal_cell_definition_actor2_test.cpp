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

double custom_value(Cell* pCell, const std::string& name)
{
    const int idx = pCell->custom_data.find_variable_index(name);
    assert(idx >= 0);
    return pCell->custom_data[idx];
}

void set_custom(Cell* pCell, const std::string& name, double value)
{
    const int idx = pCell->custom_data.find_variable_index(name);
    assert(idx >= 0);
    pCell->custom_data[idx] = value;
}

bool nearly_equal(double a, double b, double tol = 1e-12)
{
    return std::fabs(a - b) <= tol;
}

} // namespace

int main()
{
    const bool xml_status = load_PhysiCell_config_file("config/PhysiCell_settings.xml");
    assert(xml_status);

    setup_microenvironment();
    create_cell_container_for_microenvironment(microenvironment, 30.0);
    create_cell_types();

    Cell_Definition* pStroma = find_cell_definition("stromal_cell");
    assert(pStroma != NULL);
    assert(drug_index >= 0);
    assert(tgfb_index >= 0);
    assert(shh_index >= 0);

    // 1) Create PSC and check baseline activation state.
    Cell* psc = create_cell(*pStroma);
    psc->assign_position(std::vector<double>{100.0, 100.0, 0.0});
    assert(nearly_equal(custom_value(psc, "acta2_active"), 0.0));
    assert(nearly_equal(custom_value(psc, "gli1_active"), 0.0));
    assert(nearly_equal(custom_value(psc, "activation_mode"), 0.0));

    // 2) SMAD4 is WT and fixed for stroma.
    assert(nearly_equal(custom_value(psc, "gene_state_SMAD4"), 0.0));

    // 3) Death rate is exactly zero for stromal cells.
    const int apoptosis_idx =
        psc->phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);
    assert(apoptosis_idx >= 0);
    assert(nearly_equal(psc->phenotype.death.rates[apoptosis_idx], 0.0));

    // Forbidden actor-2 variables must not exist on stromal cells.
    const std::vector<std::string> forbidden = {
        "zeb1_active", "cdh1_expressed", "mmp2_active", "emt_extent",
        "nrf2_active", "abcb1_active", "intracellular_drug", "time_since_drug_exposure"
    };
    for (size_t i = 0; i < forbidden.size(); ++i)
    {
        assert(psc->custom_data.find_variable_index(forbidden[i]) < 0);
    }

    // 4) Drug-indifferent check under persistent high drug.
    std::vector<double>& psc_rho = psc->nearest_density_vector();
    psc_rho[drug_index] = 1.0;
    for (int step = 0; step < 100; ++step)
    {
        custom_function(psc, psc->phenotype, 1.0);
        assert(!psc->phenotype.death.dead);
        assert(nearly_equal(psc->phenotype.death.rates[apoptosis_idx], 0.0));
        assert(nearly_equal(psc->phenotype.secretion.uptake_rates[drug_index], 0.0));
    }
    assert(psc->custom_data.find_variable_index("intracellular_drug") < 0);

    // 5) Activate PSC -> CAF.
    parameters.doubles("tgfb_activation_threshold") = 0.3;
    parameters.doubles("shh_activation_threshold") = 0.05;
    psc_rho[tgfb_index] = 0.5;
    psc_rho[shh_index] = 0.3;
    module3_stromal_activation(psc, psc->phenotype, 1.0, ModulePhase::SENSING);
    assert(nearly_equal(custom_value(psc, "acta2_active"), 1.0));
    assert(nearly_equal(custom_value(psc, "activation_mode"), 1.0));

    // 6) Attempt deactivation by removing activating signals; must stay ON.
    psc_rho[tgfb_index] = 0.0;
    psc_rho[shh_index] = 0.0;
    module3_stromal_activation(psc, psc->phenotype, 1.0, ModulePhase::SENSING);
    assert(nearly_equal(custom_value(psc, "acta2_active"), 1.0));
    assert(nearly_equal(custom_value(psc, "activation_mode"), 1.0));

    // 7) CAF in drug field: death remains zero, TGFB active, SHH always zero.
    Cell* caf = create_cell(*pStroma);
    caf->assign_position(std::vector<double>{140.0, 100.0, 0.0});
    std::vector<double>& caf_rho = caf->nearest_density_vector();
    caf_rho[tgfb_index] = 0.6;
    caf_rho[shh_index] = 0.2;
    caf_rho[drug_index] = 1.0;
    module3_stromal_activation(caf, caf->phenotype, 1.0, ModulePhase::SENSING);
    module4_proliferation_death(caf, caf->phenotype, 1.0, ModulePhase::DECISION);
    module2_paracrine_secretion(caf, caf->phenotype, 1.0, ModulePhase::WRITE);

    const int caf_apoptosis_idx =
        caf->phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);
    assert(caf_apoptosis_idx >= 0);
    assert(nearly_equal(caf->phenotype.death.rates[caf_apoptosis_idx], 0.0));
    assert(nearly_equal(custom_value(caf, "tgfb_secretion_active"), 1.0));
    assert(caf->phenotype.secretion.secretion_rates[tgfb_index] > 0.0);
    assert(nearly_equal(caf->phenotype.secretion.secretion_rates[shh_index], 0.0));

    // 8) Rule 31 — CAFs are drug-indifferent (10 CAFs, high drug, 200 steps).
    parameters.doubles("caf_proliferation_rate") = 0.0;
    std::vector<Cell*> caf_batch_drug;
    caf_batch_drug.reserve(10);
    for (int i = 0; i < 10; ++i)
    {
        Cell* p = create_cell(*pStroma);
        p->assign_position(std::vector<double>{200.0 + 5.0 * i, 120.0, 0.0});
        set_custom(p, "acta2_active", 1.0);
        set_custom(p, "activation_mode", 1.0);
        set_custom(p, "gli1_active", 0.0);
        set_custom(p, "tgfb_secretion_active", 1.0);
        p->nearest_density_vector()[drug_index] = 1.0;
        caf_batch_drug.push_back(p);
    }

    for (int step = 0; step < 200; ++step)
    {
        for (size_t i = 0; i < caf_batch_drug.size(); ++i)
        {
            Cell* p = caf_batch_drug[i];
            custom_function(p, p->phenotype, 1.0);
            assert(!p->phenotype.death.dead);
            assert(nearly_equal(p->phenotype.death.rates[apoptosis_idx], 0.0));
            assert(nearly_equal(p->phenotype.secretion.uptake_rates[drug_index], 0.0));
        }
    }

    int alive_drug_batch = 0;
    for (size_t i = 0; i < caf_batch_drug.size(); ++i)
    {
        if (!caf_batch_drug[i]->phenotype.death.dead) ++alive_drug_batch;
    }
    assert(alive_drug_batch == 10);
    std::cout << "PASS Rule31_caf_drug_indifferent" << std::endl;

    // 9) Rule 32 — CAF death rate remains zero over long horizon (20 CAFs, 1000 steps).
    parameters.doubles("caf_proliferation_rate") = 0.0;
    std::vector<Cell*> caf_batch_long;
    caf_batch_long.reserve(20);
    for (int i = 0; i < 20; ++i)
    {
        Cell* p = create_cell(*pStroma);
        p->assign_position(std::vector<double>{250.0 + 5.0 * i, 160.0, 0.0});
        set_custom(p, "acta2_active", 1.0);
        set_custom(p, "activation_mode", 1.0);
        set_custom(p, "gli1_active", 0.0);
        set_custom(p, "tgfb_secretion_active", 1.0);
        p->nearest_density_vector()[drug_index] = 0.0;
        caf_batch_long.push_back(p);
    }

    for (int step = 0; step < 1000; ++step)
    {
        for (size_t i = 0; i < caf_batch_long.size(); ++i)
        {
            Cell* p = caf_batch_long[i];
            custom_function(p, p->phenotype, 1.0);
            assert(!p->phenotype.death.dead);
            assert(nearly_equal(p->phenotype.death.rates[apoptosis_idx], 0.0));
        }
    }

    int alive_long_batch = 0;
    for (size_t i = 0; i < caf_batch_long.size(); ++i)
    {
        if (!caf_batch_long[i]->phenotype.death.dead) ++alive_long_batch;
    }
    assert(alive_long_batch == 20);
    std::cout << "PASS Rule32_caf_zero_death" << std::endl;

    std::cout << "PASS stromal_cell_definition_actor2_test" << std::endl;
    return 0;
}
