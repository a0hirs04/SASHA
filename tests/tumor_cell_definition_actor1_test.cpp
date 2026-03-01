#include <cassert>
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

} // namespace

int main()
{
    const bool xml_status = load_PhysiCell_config_file("config/PhysiCell_settings.xml");
    assert(xml_status);

    setup_microenvironment();
    create_cell_container_for_microenvironment(microenvironment, 30.0);
    create_cell_types();

    Cell_Definition* pTumorDef = find_cell_definition("tumor_cell");
    assert(pTumorDef != NULL);

    Cell* parent = create_cell(*pTumorDef);
    assert(parent != NULL);

    // Explicitly stamp founder identity for direct-create test setup.
    set_custom(parent, "cell_id", static_cast<double>(parent->ID));
    set_custom(parent, "parent_id", -1.0);

    const std::vector<std::pair<std::string, double> > canonical_gene_states = {
        {"gene_state_KRAS", 1.0},
        {"gene_state_TP53", 2.0},
        {"gene_state_CDKN2A", 2.0},
        {"gene_state_SMAD4", 2.0},
        {"gene_state_BCL_XL", 1.0},
        {"gene_state_MYC", 0.0},
        {"gene_state_EGFR", 0.0},
        {"gene_state_RB1", 0.0},
        {"gene_state_BAX", 0.0},
        {"gene_state_ZEB1", 0.0},
        {"gene_state_CDH1", 0.0},
        {"gene_state_MMP2", 0.0},
        {"gene_state_TGFB1", 0.0},
        {"gene_state_SHH", 0.0},
        {"gene_state_HIF1A", 0.0},
        {"gene_state_NRF2", 0.0},
        {"gene_state_ABCB1", 0.0},
        {"gene_state_HAS2", 0.0},
        {"gene_state_COL1A1", 0.0},
        {"gene_state_GLI1", 0.0},
        {"gene_state_ACTA2", 0.0}
    };
    assert(canonical_gene_states.size() == 21);

    for (size_t i = 0; i < canonical_gene_states.size(); ++i)
    {
        assert(custom_value(parent, canonical_gene_states[i].first) == canonical_gene_states[i].second);
    }

    // Required Actor 1 defaults.
    assert(custom_value(parent, "hif1a_active") == 0.0);
    assert(custom_value(parent, "zeb1_active") == 0.0);
    assert(custom_value(parent, "cdh1_expressed") == 1.0);
    assert(custom_value(parent, "nrf2_active") == 0.0);
    assert(custom_value(parent, "abcb1_active") == 0.0);
    assert(custom_value(parent, "mmp2_active") == 0.0);
    assert(custom_value(parent, "intracellular_drug") == 0.0);
    assert(custom_value(parent, "time_since_drug_exposure") == -1.0);

    // Tumor cells must not register stromal-only activity flags.
    assert(parent->custom_data.find_variable_index("acta2_active") < 0);
    assert(parent->custom_data.find_variable_index("gli1_active") < 0);

    // Put parent in a non-default state; daughter must reset by division callback.
    set_custom(parent, "hif1a_active", 1.0);
    set_custom(parent, "zeb1_active", 1.0);
    set_custom(parent, "cdh1_expressed", 0.0);
    set_custom(parent, "nrf2_active", 1.0);
    set_custom(parent, "abcb1_active", 1.0);
    set_custom(parent, "mmp2_active", 1.0);
    set_custom(parent, "intracellular_drug", 0.72);
    set_custom(parent, "time_since_drug_exposure", 18.0);
    set_custom(parent, "emt_extent", 3.0);
    set_custom(parent, "time_alive", 120.0);
    set_custom(parent, "mechanical_pressure", 2.0);

    const double parent_cell_id = custom_value(parent, "cell_id");
    assert(parent_cell_id == static_cast<double>(parent->ID));

    Cell* daughter = parent->divide();
    assert(daughter != NULL);

    for (size_t i = 0; i < canonical_gene_states.size(); ++i)
    {
        const double parent_gene = custom_value(parent, canonical_gene_states[i].first);
        const double daughter_gene = custom_value(daughter, canonical_gene_states[i].first);
        assert(daughter_gene == parent_gene);
    }

    assert(custom_value(daughter, "hif1a_active") == 0.0);
    assert(custom_value(daughter, "zeb1_active") == 0.0);
    assert(custom_value(daughter, "cdh1_expressed") == 1.0);
    assert(custom_value(daughter, "nrf2_active") == 0.0);
    assert(custom_value(daughter, "abcb1_active") == 0.0);
    assert(custom_value(daughter, "mmp2_active") == 0.0);
    assert(custom_value(daughter, "intracellular_drug") == 0.0);
    assert(custom_value(daughter, "time_since_drug_exposure") == -1.0);
    assert(custom_value(daughter, "emt_extent") == 0.0);
    assert(custom_value(daughter, "time_alive") == 0.0);

    const double daughter_cell_id = custom_value(daughter, "cell_id");
    const double daughter_parent_id = custom_value(daughter, "parent_id");
    assert(daughter_cell_id == static_cast<double>(daughter->ID));
    assert(daughter_cell_id != parent_cell_id);
    assert(daughter_parent_id == parent_cell_id);

    std::cout << "PASS tumor_cell_definition_actor1_test" << std::endl;
    return 0;
}
