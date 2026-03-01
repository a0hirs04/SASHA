#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "feedback_loop_common.h"

using namespace BioFVM;
using namespace PhysiCell;
using namespace feedback_loop_common;

namespace
{

bool nearly_equal(double a, double b, double eps = 1e-10)
{
    return std::fabs(a - b) <= eps;
}

int idx(Cell* pCell, const std::string& name)
{
    const int i = pCell->custom_data.find_variable_index(name);
    assert(i >= 0);
    return i;
}

double cget(Cell* pCell, const std::string& name)
{
    return pCell->custom_data[idx(pCell, name)];
}

void cset_if_present(Cell* pCell, const std::string& name, double value)
{
    const int i = pCell->custom_data.find_variable_index(name);
    if (i >= 0) pCell->custom_data[i] = value;
}

void clear_all_cells()
{
    while (!all_cells->empty())
    {
        delete_cell((*all_cells)[all_cells->size() - 1]);
    }
}

void set_uniform_field(int substrate_index, double value)
{
    const int n_vox = static_cast<int>(microenvironment.number_of_voxels());
    for (int n = 0; n < n_vox; ++n)
    {
        std::vector<double>& rho = microenvironment.density_vector(n);
        if (substrate_index >= 0 && substrate_index < static_cast<int>(rho.size()))
        {
            rho[substrate_index] = value;
        }
    }
}

void set_local_field(Cell* pCell, int substrate_index, double value)
{
    std::vector<double>& rho = pCell->nearest_density_vector();
    if (substrate_index >= 0 && substrate_index < static_cast<int>(rho.size()))
    {
        rho[substrate_index] = value;
    }
}

void neutralize_tumor_state(Cell* pCell)
{
    cset_if_present(pCell, "hif1a_active", 0.0);
    cset_if_present(pCell, "zeb1_active", 0.0);
    cset_if_present(pCell, "abcb1_active", 0.0);
    cset_if_present(pCell, "nrf2_active", 0.0);
    cset_if_present(pCell, "mmp2_active", 0.0);
    cset_if_present(pCell, "intracellular_drug", 0.0);
    cset_if_present(pCell, "time_since_drug_exposure", -1.0);
    cset_if_present(pCell, "mechanical_pressure", 0.0);
    cset_if_present(pCell, "SMAD4", 0.0);
    set_local_field(pCell, tgfb_index, 0.0);
    set_local_field(pCell, shh_index, 0.0);
    set_local_field(pCell, drug_index, 0.0);
    set_local_field(pCell, oxygen_index, 0.08);
}

int live_tumor_count(int tumor_type)
{
    int n = 0;
    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->phenotype.death.dead) continue;
        if (pCell->type == tumor_type) ++n;
    }
    return n;
}

int live_stromal_count(int stromal_type)
{
    int n = 0;
    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->phenotype.death.dead) continue;
        if (pCell->type == stromal_type) ++n;
    }
    return n;
}

double mean_live_custom(int type_id, const std::string& var_name)
{
    double sum = 0.0;
    int n = 0;
    for (size_t i = 0; i < all_cells->size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        if (pCell == NULL || pCell->phenotype.death.dead) continue;
        if (pCell->type != type_id) continue;
        const int vi = pCell->custom_data.find_variable_index(var_name);
        if (vi < 0) continue;
        sum += pCell->custom_data[vi];
        ++n;
    }
    if (n == 0) return 0.0;
    return sum / static_cast<double>(n);
}

void run_population_steps(double& t, int steps, double dt)
{
    Cell_Container* cc = static_cast<Cell_Container*>(microenvironment.agent_container);
    for (int k = 0; k < steps; ++k)
    {
        microenvironment.simulate_diffusion_decay(dt);
        cc->update_all_cells(t, dt, dt, dt);
        t += dt;
    }
}

void print_pass(const std::string& name)
{
    std::cout << "PASS " << name << std::endl;
}

} // namespace

int main()
{
    initialize_world();
    SeedRandom(0);

    assert(oxygen_index >= 0);
    assert(tgfb_index >= 0);
    assert(shh_index >= 0);
    assert(drug_index >= 0);
    assert(ecm_index >= 0);

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    Cell_Definition* pStroma = find_cell_definition("stromal_cell");
    assert(pTumor != NULL);
    assert(pStroma != NULL);

    // Shared defaults.
    parameters.doubles("base_proliferation_rate") = 0.01;
    parameters.doubles("apoptosis_resistance") = 0.85;
    parameters.doubles("go_grow_penalty") = 0.3;
    parameters.doubles("hypoxia_proliferation_modifier") = 0.2;
    parameters.doubles("hypoxia_death_resistance_bonus") = 0.1;
    parameters.doubles("contact_inhibition_threshold") = 5.0;
    parameters.doubles("drug_kill_coefficient") = 0.01;
    parameters.doubles("efflux_drug_reduction") = 0.7;
    parameters.doubles("tgfb_brake_sensitivity") = 0.0;

    parameters.doubles("drug_uptake_rate") = 0.01;
    parameters.doubles("drug_stress_threshold") = 0.03;
    parameters.doubles("efflux_induction_delay") = 12.0;
    parameters.doubles("efflux_strength") = 0.1;
    parameters.doubles("nrf2_decay_rate") = 0.01;
    parameters.doubles("drug_natural_decay_rate") = 0.005;
    parameters.doubles("hif1a_nrf2_priming_bonus") = 0.02;

    parameters.doubles("tgfb_secretion_rate") = 0.5;
    parameters.doubles("shh_secretion_rate") = 0.3;
    parameters.doubles("hif1a_tgfb_amplification_factor") = 1.5;

    parameters.doubles("emt_induction_threshold") = 0.4;
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

    // ---------------------------------------------------------------------
    // Rule 1: Genotype determines division rate + base rate matches Knob 3a.
    // ---------------------------------------------------------------------
    {
        clear_all_cells();
        Cell* fast = create_cell(*pTumor);
        Cell* slow = create_cell(*pTumor);
        Cell* canonical = create_cell(*pTumor);

        fast->assign_position(std::vector<double>{-100.0, 0.0, 0.0});
        slow->assign_position(std::vector<double>{0.0, 0.0, 0.0});
        canonical->assign_position(std::vector<double>{100.0, 0.0, 0.0});

        neutralize_tumor_state(fast);
        neutralize_tumor_state(slow);
        neutralize_tumor_state(canonical);

        cset_if_present(fast, "gene_state_KRAS", 1.0);
        cset_if_present(fast, "gene_state_MYC", 1.0);
        cset_if_present(slow, "gene_state_KRAS", 0.0);
        cset_if_present(slow, "gene_state_MYC", 0.0);
        cset_if_present(canonical, "gene_state_KRAS", 1.0);
        cset_if_present(canonical, "gene_state_MYC", 0.0);

        module4_proliferation_death(fast, fast->phenotype, 6.0, ModulePhase::DECISION);
        module4_proliferation_death(slow, slow->phenotype, 6.0, ModulePhase::DECISION);
        module4_proliferation_death(canonical, canonical->phenotype, 6.0, ModulePhase::DECISION);

        const double p_fast = fast->phenotype.cycle.data.transition_rate(0, 0);
        const double p_slow = slow->phenotype.cycle.data.transition_rate(0, 0);
        const double p_canon = canonical->phenotype.cycle.data.transition_rate(0, 0);

        assert(p_fast > p_slow);
        assert(nearly_equal(p_canon, parameters.doubles("base_proliferation_rate")));
        print_pass("Rule 1");
    }

    // ---------------------------------------------------------------------
    // Rule 2: Mechanical compression halts division, release resumes.
    // ---------------------------------------------------------------------
    {
        clear_all_cells();
        Cell* p = create_cell(*pTumor);
        p->assign_position(std::vector<double>{0.0, 0.0, 0.0});
        neutralize_tumor_state(p);
        cset_if_present(p, "mechanical_pressure", 6.0);
        module4_proliferation_death(p, p->phenotype, 6.0, ModulePhase::DECISION);
        assert(nearly_equal(p->phenotype.cycle.data.transition_rate(0, 0), 0.0));

        cset_if_present(p, "mechanical_pressure", 0.0);
        module4_proliferation_death(p, p->phenotype, 6.0, ModulePhase::DECISION);
        assert(p->phenotype.cycle.data.transition_rate(0, 0) > 0.0);
        print_pass("Rule 2");
    }

    // ---------------------------------------------------------------------
    // Rule 3: EMT (ZEB1 ON) slows division via go-grow penalty.
    // ---------------------------------------------------------------------
    {
        clear_all_cells();
        Cell* off = create_cell(*pTumor);
        Cell* on = create_cell(*pTumor);
        off->assign_position(std::vector<double>{-20.0, 0.0, 0.0});
        on->assign_position(std::vector<double>{20.0, 0.0, 0.0});
        neutralize_tumor_state(off);
        neutralize_tumor_state(on);

        cset_if_present(on, "zeb1_active", 1.0);
        module4_proliferation_death(off, off->phenotype, 6.0, ModulePhase::DECISION);
        module4_proliferation_death(on, on->phenotype, 6.0, ModulePhase::DECISION);

        const double p_off = off->phenotype.cycle.data.transition_rate(0, 0);
        const double p_on = on->phenotype.cycle.data.transition_rate(0, 0);
        assert(p_on < p_off);
        assert(nearly_equal(p_on, p_off * (1.0 - parameters.doubles("go_grow_penalty"))));
        print_pass("Rule 3");
    }

    // ---------------------------------------------------------------------
    // Rule 4 + Rule 5: Division inheritance/reset behavior.
    // ---------------------------------------------------------------------
    {
        clear_all_cells();
        Cell* parent = create_cell(*pTumor);
        parent->assign_position(std::vector<double>{0.0, 0.0, 0.0});
        cset_if_present(parent, "cell_id", static_cast<double>(parent->ID));
        cset_if_present(parent, "parent_id", -1.0);

        const std::vector<std::string> gene_names = {
            "gene_state_KRAS", "gene_state_TP53", "gene_state_CDKN2A", "gene_state_SMAD4",
            "gene_state_BCL_XL", "gene_state_MYC", "gene_state_EGFR", "gene_state_RB1",
            "gene_state_BAX", "gene_state_ZEB1", "gene_state_CDH1", "gene_state_MMP2",
            "gene_state_TGFB1", "gene_state_SHH", "gene_state_HIF1A", "gene_state_NRF2",
            "gene_state_ABCB1", "gene_state_HAS2", "gene_state_COL1A1", "gene_state_GLI1",
            "gene_state_ACTA2"
        };
        assert(gene_names.size() == 21);

        for (size_t i = 0; i < gene_names.size(); ++i)
        {
            cset_if_present(parent, gene_names[i], (i % 3 == 0) ? 2.0 : ((i % 3 == 1) ? 1.0 : 0.0));
        }

        cset_if_present(parent, "hif1a_active", 1.0);
        cset_if_present(parent, "zeb1_active", 1.0);
        cset_if_present(parent, "cdh1_expressed", 0.0);
        cset_if_present(parent, "nrf2_active", 1.0);
        cset_if_present(parent, "abcb1_active", 1.0);
        cset_if_present(parent, "mmp2_active", 1.0);
        cset_if_present(parent, "intracellular_drug", 0.6);
        cset_if_present(parent, "time_since_drug_exposure", 50.0);

        const double parent_cell_id = static_cast<double>(parent->ID);
        Cell* daughter = parent->divide();
        assert(daughter != NULL);

        for (size_t i = 0; i < gene_names.size(); ++i)
        {
            assert(nearly_equal(cget(daughter, gene_names[i]), cget(parent, gene_names[i])));
        }

        assert(nearly_equal(cget(daughter, "hif1a_active"), 0.0));
        assert(nearly_equal(cget(daughter, "zeb1_active"), 0.0));
        assert(nearly_equal(cget(daughter, "cdh1_expressed"), 1.0));
        assert(nearly_equal(cget(daughter, "nrf2_active"), 0.0));
        assert(nearly_equal(cget(daughter, "abcb1_active"), 0.0));
        assert(nearly_equal(cget(daughter, "mmp2_active"), 0.0));

        assert(nearly_equal(cget(daughter, "intracellular_drug"), 0.0));
        assert(nearly_equal(cget(daughter, "time_since_drug_exposure"), -1.0));

        const double callback_parent_cell_id = cget(parent, "cell_id");
        assert(nearly_equal(callback_parent_cell_id, parent_cell_id));
        assert(!nearly_equal(cget(daughter, "cell_id"), callback_parent_cell_id));
        assert(nearly_equal(cget(daughter, "parent_id"), callback_parent_cell_id));

        print_pass("Rule 4");
        print_pass("Rule 5");
    }

    // ---------------------------------------------------------------------
    // Rule 6: Baseline apoptosis low, drug increases apoptosis.
    // ---------------------------------------------------------------------
    {
        clear_all_cells();
        Cell* p = create_cell(*pTumor);
        p->assign_position(std::vector<double>{0.0, 0.0, 0.0});
        neutralize_tumor_state(p);

        const int apop_idx = p->phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);
        assert(apop_idx >= 0);

        module4_proliferation_death(p, p->phenotype, 6.0, ModulePhase::DECISION);
        const double baseline = p->phenotype.death.rates[apop_idx];

        cset_if_present(p, "intracellular_drug", 0.5);
        module4_proliferation_death(p, p->phenotype, 6.0, ModulePhase::DECISION);
        const double with_drug = p->phenotype.death.rates[apop_idx];

        assert(with_drug > baseline);
        print_pass("Rule 6");
    }

    // ---------------------------------------------------------------------
    // Rule 7: Dead cells removed; no secretion/uptake/pressure contribution.
    // ---------------------------------------------------------------------
    {
        clear_all_cells();

        Cell* dying = create_cell(*pTumor);
        const std::vector<double> pos = {15.0, 15.0, 0.0};
        dying->assign_position(pos);
        neutralize_tumor_state(dying);
        cset_if_present(dying, "mechanical_pressure", 2.5);
        module2_paracrine_secretion(dying, dying->phenotype, 6.0, ModulePhase::WRITE);

        const int apop_idx = dying->phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);
        assert(apop_idx >= 0);
        dying->start_death(apop_idx);
        custom_function(dying, dying->phenotype, 6.0);

        assert(nearly_equal(dying->phenotype.secretion.uptake_rates[oxygen_index], 0.0));
        assert(nearly_equal(dying->phenotype.secretion.secretion_rates[tgfb_index], 0.0));
        assert(nearly_equal(dying->phenotype.secretion.secretion_rates[shh_index], 0.0));
        assert(nearly_equal(cget(dying, "mechanical_pressure"), 0.0));

        const size_t before = all_cells->size();
        dying->die();
        assert(all_cells->size() + 1 == before);

        // Space is available: place a new cell at same position.
        Cell* replacement = create_cell(*pTumor);
        const bool assigned = replacement->assign_position(pos);
        assert(assigned);
        assert(!replacement->phenotype.death.dead);

        print_pass("Rule 7");
    }

    // ---------------------------------------------------------------------
    // Rule 8: First drug pulse is more dangerous than second pulse.
    // ---------------------------------------------------------------------
    {
        clear_all_cells();
        parameters.doubles("base_proliferation_rate") = 0.0; // avoid births
        parameters.doubles("apoptosis_resistance") = 0.98;
        parameters.doubles("drug_kill_coefficient") = 0.0020;
        parameters.doubles("drug_uptake_rate") = 0.015;
        parameters.doubles("efflux_strength") = 0.2;
        parameters.doubles("efflux_induction_delay") = 240.0;
        parameters.doubles("drug_natural_decay_rate") = 0.004;
        parameters.doubles("drug_stress_threshold") = 0.03;

        for (int i = 0; i < 20; ++i)
        {
            Cell* p = create_cell(*pTumor);
            const double x = -60.0 + 6.0 * static_cast<double>(i);
            p->assign_position(std::vector<double>{x, 0.0, 0.0});
            neutralize_tumor_state(p);
        }

        Cell_Container* cc = static_cast<Cell_Container*>(microenvironment.agent_container);
        double t = 0.0;
        const double dt = 6.0;

        const int tumor_type = pTumor->type;
        const int initial = live_tumor_count(tumor_type);

        // Pulse 1: 100 steps with drug.
        for (int step = 0; step < 100; ++step)
        {
            for (size_t i = 0; i < all_cells->size(); ++i)
            {
                Cell* p = (*all_cells)[i];
                if (p == NULL || p->phenotype.death.dead || p->type != tumor_type) continue;
                set_local_field(p, drug_index, 1.0);
            }
            cc->update_all_cells(t, dt, dt, dt);
            t += dt;
        }
        const int survivors_after_pulse1 = live_tumor_count(tumor_type);
        assert(survivors_after_pulse1 > 0);

        // Washout: 50 steps without drug.
        for (int step = 0; step < 50; ++step)
        {
            for (size_t i = 0; i < all_cells->size(); ++i)
            {
                Cell* p = (*all_cells)[i];
                if (p == NULL || p->phenotype.death.dead || p->type != tumor_type) continue;
                set_local_field(p, drug_index, 0.0);
            }
            cc->update_all_cells(t, dt, dt, dt);
            t += dt;
        }

        // Pulse 2: 100 steps with drug.
        for (int step = 0; step < 100; ++step)
        {
            for (size_t i = 0; i < all_cells->size(); ++i)
            {
                Cell* p = (*all_cells)[i];
                if (p == NULL || p->phenotype.death.dead || p->type != tumor_type) continue;
                set_local_field(p, drug_index, 1.0);
            }
            cc->update_all_cells(t, dt, dt, dt);
            t += dt;
        }
        const int survivors_after_pulse2 = live_tumor_count(tumor_type);

        const int deaths_pulse1 = initial - survivors_after_pulse1;
        const int deaths_pulse2 = survivors_after_pulse1 - survivors_after_pulse2;

        const double abcb1_after_pulse2 = mean_live_custom(tumor_type, "abcb1_active");

        std::cerr << "Rule8 metrics"
                  << " initial=" << initial
                  << " after_pulse1=" << survivors_after_pulse1
                  << " after_pulse2=" << survivors_after_pulse2
                  << " deaths_pulse1=" << deaths_pulse1
                  << " deaths_pulse2=" << deaths_pulse2
                  << " mean_abcb1_survivors=" << abcb1_after_pulse2
                  << std::endl;

        const bool pass_first_more_dangerous = (deaths_pulse1 > deaths_pulse2);
        const bool pass_memory_present = (abcb1_after_pulse2 > 0.0);
        assert(pass_first_more_dangerous);
        assert(pass_memory_present);
        print_pass("Rule 8");

        // restore default for subsequent tests
        parameters.doubles("base_proliferation_rate") = 0.01;
        parameters.doubles("apoptosis_resistance") = 0.85;
        parameters.doubles("drug_kill_coefficient") = 0.01;
    }

    // ---------------------------------------------------------------------
    // Rule 9: Gene states never change over time.
    // ---------------------------------------------------------------------
    {
        clear_all_cells();
        parameters.doubles("base_proliferation_rate") = 0.0;
        parameters.doubles("apoptosis_resistance") = 1.0;
        parameters.doubles("drug_kill_coefficient") = 0.0;

        const std::vector<std::string> gene_names = {
            "gene_state_KRAS", "gene_state_TP53", "gene_state_CDKN2A", "gene_state_SMAD4",
            "gene_state_BCL_XL", "gene_state_MYC", "gene_state_EGFR", "gene_state_RB1",
            "gene_state_BAX", "gene_state_ZEB1", "gene_state_CDH1", "gene_state_MMP2",
            "gene_state_TGFB1", "gene_state_SHH", "gene_state_HIF1A", "gene_state_NRF2",
            "gene_state_ABCB1", "gene_state_HAS2", "gene_state_COL1A1", "gene_state_GLI1",
            "gene_state_ACTA2"
        };

        std::vector<Cell*> cells;
        for (int i = 0; i < 5; ++i)
        {
            Cell* p = create_cell(*pTumor);
            p->assign_position(std::vector<double>{static_cast<double>(i) * 30.0, 0.0, 0.0});
            neutralize_tumor_state(p);
            cells.push_back(p);
        }

        std::map<int, std::vector<double> > baseline;
        for (size_t i = 0; i < cells.size(); ++i)
        {
            std::vector<double> vals;
            for (size_t g = 0; g < gene_names.size(); ++g)
            {
                vals.push_back(cget(cells[i], gene_names[g]));
            }
            baseline[cells[i]->ID] = vals;
        }

        for (int step = 0; step < 500; ++step)
        {
            for (size_t i = 0; i < cells.size(); ++i)
            {
                custom_function(cells[i], cells[i]->phenotype, 1.0);
            }
        }

        for (size_t i = 0; i < cells.size(); ++i)
        {
            Cell* p = cells[i];
            assert(!p->phenotype.death.dead);
            const std::vector<double>& base = baseline[p->ID];
            for (size_t g = 0; g < gene_names.size(); ++g)
            {
                assert(nearly_equal(cget(p, gene_names[g]), base[g]));
            }
        }

        print_pass("Rule 9");

        parameters.doubles("base_proliferation_rate") = 0.01;
        parameters.doubles("apoptosis_resistance") = 0.85;
        parameters.doubles("drug_kill_coefficient") = 0.01;
    }

    // ---------------------------------------------------------------------
    // Rule 10: HIF1A adaptation is reversible.
    // ---------------------------------------------------------------------
    {
        clear_all_cells();
        Cell* p = create_cell(*pTumor);
        p->assign_position(std::vector<double>{0.0, 0.0, 0.0});

        set_local_field(p, oxygen_index, 0.01);
        module1_oxygen_sensing(p, p->phenotype, 6.0, ModulePhase::SENSING);
        assert(nearly_equal(cget(p, "hif1a_active"), 1.0));

        set_local_field(p, oxygen_index, 0.08);
        for (int i = 0; i < 10; ++i)
        {
            module1_oxygen_sensing(p, p->phenotype, 6.0, ModulePhase::SENSING);
        }
        assert(nearly_equal(cget(p, "hif1a_active"), 0.0));
        print_pass("Rule 10");
    }

    // ---------------------------------------------------------------------
    // Rule 11: ABCB1 memory persists while NRF2 decays under no drug.
    // ---------------------------------------------------------------------
    {
        clear_all_cells();
        Cell* p = create_cell(*pTumor);
        p->assign_position(std::vector<double>{0.0, 0.0, 0.0});
        neutralize_tumor_state(p);

        parameters.doubles("efflux_induction_delay") = 6.0;
        parameters.doubles("nrf2_decay_rate") = 0.01;

        // Activate with drug for several steps.
        for (int i = 0; i < 8; ++i)
        {
            set_local_field(p, drug_index, 1.0);
            module7_drug_response(p, p->phenotype, 6.0, ModulePhase::SENSING);
        }
        assert(cget(p, "abcb1_active") >= 1.0);

        // Remove drug and allow decay.
        for (int i = 0; i < 50; ++i)
        {
            set_local_field(p, drug_index, 0.0);
            module7_drug_response(p, p->phenotype, 6.0, ModulePhase::SENSING);
        }

        assert(cget(p, "nrf2_active") < 1.0);
        assert(cget(p, "abcb1_active") >= 1.0);
        print_pass("Rule 11");
    }

    // ---------------------------------------------------------------------
    // Rule 13 + Rule 14: Constitutive secretion + HIF1A amplifies TGFb only.
    // ---------------------------------------------------------------------
    {
        clear_all_cells();
        Cell* p = create_cell(*pTumor);
        p->assign_position(std::vector<double>{0.0, 0.0, 0.0});
        neutralize_tumor_state(p);

        // Stress conditions: drug, hypoxia, compression.
        set_local_field(p, drug_index, 1.0);
        set_local_field(p, oxygen_index, 0.01);
        cset_if_present(p, "mechanical_pressure", 9.0);
        module1_oxygen_sensing(p, p->phenotype, 6.0, ModulePhase::SENSING);
        module2_paracrine_secretion(p, p->phenotype, 6.0, ModulePhase::WRITE);

        const double tgfb_stress = p->phenotype.secretion.secretion_rates[tgfb_index];
        const double shh_stress = p->phenotype.secretion.secretion_rates[shh_index];
        assert(tgfb_stress > 0.0);
        assert(shh_stress > 0.0);
        print_pass("Rule 13");

        cset_if_present(p, "hif1a_active", 0.0);
        module2_paracrine_secretion(p, p->phenotype, 6.0, ModulePhase::WRITE);
        const double tgfb_hif0 = p->phenotype.secretion.secretion_rates[tgfb_index];
        const double shh_hif0 = p->phenotype.secretion.secretion_rates[shh_index];

        cset_if_present(p, "hif1a_active", 1.0);
        module2_paracrine_secretion(p, p->phenotype, 6.0, ModulePhase::WRITE);
        const double tgfb_hif1 = p->phenotype.secretion.secretion_rates[tgfb_index];
        const double shh_hif1 = p->phenotype.secretion.secretion_rates[shh_index];

        assert(tgfb_hif1 > tgfb_hif0);
        assert(nearly_equal(shh_hif1, shh_hif0));
        print_pass("Rule 14");
    }

    // ---------------------------------------------------------------------
    // Rule 15: ECM slows diffusion (effective coefficient drop).
    // ---------------------------------------------------------------------
    {
        clear_all_cells();
        const int voxel = voxel_index_for_position(0.0, 0.0, 0.0);
        assert(voxel >= 0);

        microenvironment.density_vector(voxel)[ecm_index] = 0.0;
        set_ecm_ha_fraction(voxel, 0.8);
        update_ecm_effective_diffusion_coefficients(microenvironment);
        const double d_tgfb_no_ecm = get_effective_diffusion_coefficient(tgfb_index, voxel);

        microenvironment.density_vector(voxel)[ecm_index] = 0.8;
        set_ecm_ha_fraction(voxel, 0.8);
        update_ecm_effective_diffusion_coefficients(microenvironment);
        const double d_tgfb_with_ecm = get_effective_diffusion_coefficient(tgfb_index, voxel);

        assert(d_tgfb_with_ecm < d_tgfb_no_ecm);
        print_pass("Rule 15");
    }

    // ---------------------------------------------------------------------
    // Rule 16 + Rule 17 + Rule 18: EMT trigger, reversibility, phenotype.
    // ---------------------------------------------------------------------
    {
        clear_all_cells();
        Cell* p = create_cell(*pTumor);
        p->assign_position(std::vector<double>{0.0, 0.0, 0.0});
        neutralize_tumor_state(p);

        // E01/E08 style trigger: TGFb below threshold without HIF1A, above with HIF1A boost.
        set_local_field(p, tgfb_index, 0.3);
        cset_if_present(p, "hif1a_active", 0.0);
        module5_emt_engine(p, p->phenotype, 6.0, ModulePhase::DECISION);
        assert(nearly_equal(cget(p, "zeb1_active"), 0.0));

        cset_if_present(p, "hif1a_active", 1.0);
        module5_emt_engine(p, p->phenotype, 6.0, ModulePhase::DECISION);
        assert(nearly_equal(cget(p, "zeb1_active"), 1.0));
        assert(nearly_equal(cget(p, "cdh1_expressed"), 0.0));
        assert(nearly_equal(cget(p, "mmp2_active"), 1.0));
        assert(nearly_equal(p->phenotype.motility.migration_speed, parameters.doubles("motility_mesenchymal_med")));
        assert(nearly_equal(p->phenotype.mechanics.cell_cell_adhesion_strength, parameters.doubles("adhesion_mesenchymal_med")));

        // Reversibility
        set_local_field(p, tgfb_index, 0.0);
        cset_if_present(p, "hif1a_active", 0.0);
        module5_emt_engine(p, p->phenotype, 6.0, ModulePhase::DECISION);
        assert(nearly_equal(cget(p, "zeb1_active"), 0.0));
        assert(nearly_equal(cget(p, "cdh1_expressed"), 1.0));
        assert(nearly_equal(cget(p, "mmp2_active"), 0.0));
        assert(nearly_equal(p->phenotype.motility.migration_speed, parameters.doubles("motility_epithelial")));
        assert(nearly_equal(p->phenotype.mechanics.cell_cell_adhesion_strength, parameters.doubles("adhesion_epithelial")));

        print_pass("Rule 16");
        print_pass("Rule 17");
        print_pass("Rule 18");
    }

    // ---------------------------------------------------------------------
    // Rule 19: MMP2 ECM degradation is local.
    // ---------------------------------------------------------------------
    {
        clear_all_cells();
        Cell* p = create_cell(*pTumor);
        p->assign_position(std::vector<double>{0.0, 0.0, 0.0});
        neutralize_tumor_state(p);

        parameters.doubles("mmp2_degradation_rate") = 0.005;

        const int v0 = voxel_index_for_cell(p);
        const int v1 = voxel_index_for_position(40.0, 0.0, 0.0);
        assert(v0 >= 0 && v1 >= 0 && v0 != v1);

        microenvironment.density_vector(v0)[ecm_index] = 0.5;
        microenvironment.density_vector(v1)[ecm_index] = 0.5;
        set_ecm_ha_fraction(v0, 0.6);
        set_ecm_ha_fraction(v1, 0.6);

        cset_if_present(p, "mmp2_active", 1.0);
        module6_ecm_production(p, p->phenotype, 6.0, ModulePhase::WRITE);

        const double e0 = microenvironment.density_vector(v0)[ecm_index];
        const double e1 = microenvironment.density_vector(v1)[ecm_index];
        assert(e0 < 0.5);
        assert(nearly_equal(e1, 0.5));
        assert(e0 >= 0.0);

        print_pass("Rule 19");
    }

    // ---------------------------------------------------------------------
    // Rule 20: SMAD4 asymmetry (E08 always on, E09 WT-only).
    // ---------------------------------------------------------------------
    {
        clear_all_cells();

        parameters.doubles("base_proliferation_rate") = 0.01;
        parameters.doubles("go_grow_penalty") = 0.0;
        parameters.doubles("hypoxia_proliferation_modifier") = 0.0;
        parameters.doubles("contact_inhibition_threshold") = 1e9;
        parameters.doubles("tgfb_brake_sensitivity") = 0.5;
        parameters.doubles("hif1a_emt_boost") = 0.0;

        Cell* A = create_cell(*pTumor); // WT, no TGFb
        Cell* B = create_cell(*pTumor); // WT, high TGFb
        Cell* C = create_cell(*pTumor); // LOF, no TGFb
        Cell* D = create_cell(*pTumor); // LOF, high TGFb
        A->assign_position(std::vector<double>{-60.0, 0.0, 0.0});
        B->assign_position(std::vector<double>{-20.0, 0.0, 0.0});
        C->assign_position(std::vector<double>{20.0, 0.0, 0.0});
        D->assign_position(std::vector<double>{60.0, 0.0, 0.0});

        std::vector<Cell*> cells = {A, B, C, D};
        for (size_t i = 0; i < cells.size(); ++i)
        {
            neutralize_tumor_state(cells[i]);
            cset_if_present(cells[i], "hif1a_active", 0.0);
        }

        cset_if_present(A, "SMAD4", 1.0);
        cset_if_present(B, "SMAD4", 1.0);
        cset_if_present(C, "SMAD4", 0.0);
        cset_if_present(D, "SMAD4", 0.0);

        set_local_field(A, tgfb_index, 0.0);
        set_local_field(B, tgfb_index, 0.6);
        set_local_field(C, tgfb_index, 0.0);
        set_local_field(D, tgfb_index, 0.6);

        for (size_t i = 0; i < cells.size(); ++i)
        {
            module5_emt_engine(cells[i], cells[i]->phenotype, 6.0, ModulePhase::DECISION);
            module4_proliferation_death(cells[i], cells[i]->phenotype, 6.0, ModulePhase::DECISION);
        }

        const double base = parameters.doubles("base_proliferation_rate");
        const double pA = A->phenotype.cycle.data.transition_rate(0, 0);
        const double pB = B->phenotype.cycle.data.transition_rate(0, 0);
        const double pC = C->phenotype.cycle.data.transition_rate(0, 0);
        const double pD = D->phenotype.cycle.data.transition_rate(0, 0);

        assert(nearly_equal(pA, base));
        assert(pB < base);
        assert(nearly_equal(pC, base));
        assert(nearly_equal(pD, base));

        assert(nearly_equal(cget(A, "zeb1_active"), 0.0));
        assert(nearly_equal(cget(B, "zeb1_active"), 1.0));
        assert(nearly_equal(cget(C, "zeb1_active"), 0.0));
        assert(nearly_equal(cget(D, "zeb1_active"), 1.0));

        print_pass("Rule 20");
    }

    // ---------------------------------------------------------------------
    // Rule 21: Stromal SMAD4 always WT and immutable over time.
    // ---------------------------------------------------------------------
    {
        clear_all_cells();

        std::vector<Cell*> stromals;
        for (int i = 0; i < 8; ++i)
        {
            Cell* s = create_cell(*pStroma);
            s->assign_position(std::vector<double>{-80.0 + 20.0 * i, 0.0, 0.0});
            stromals.push_back(s);
        }

        for (size_t i = 0; i < stromals.size(); ++i)
        {
            assert(nearly_equal(cget(stromals[i], "gene_state_SMAD4"), 0.0));
        }

        for (int step = 0; step < 100; ++step)
        {
            // Attempt mutation; runtime should re-lock to WT.
            cset_if_present(stromals.front(), "gene_state_SMAD4", 2.0);

            for (size_t i = 0; i < stromals.size(); ++i)
            {
                set_local_field(stromals[i], tgfb_index, 0.8);
                set_local_field(stromals[i], shh_index, 0.5);
                set_local_field(stromals[i], drug_index, 1.0);
                custom_function(stromals[i], stromals[i]->phenotype, 6.0);
                assert(nearly_equal(cget(stromals[i], "gene_state_SMAD4"), 0.0));
            }
        }

        print_pass("Rule 21");
    }

    std::cout << "PASS test_tumor_rules" << std::endl;
    return 0;
}
