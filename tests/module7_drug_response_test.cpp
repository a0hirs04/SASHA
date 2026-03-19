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

int voxel_index_for_cell(Cell* pCell)
{
    std::vector<double> pos = pCell->position;
    return microenvironment.nearest_voxel_index(pos);
}

void set_local_drug(Cell* pCell, double value)
{
    const int voxel = voxel_index_for_cell(pCell);
    microenvironment.density_vector(voxel)[drug_index] = value;
}

} // namespace

int main()
{
    const bool xml_status = load_PhysiCell_config_file("config/PhysiCell_settings.xml");
    assert(xml_status);

    setup_microenvironment();
    create_cell_container_for_microenvironment(microenvironment, 30.0);
    create_cell_types();

    assert(drug_index >= 0);

    parameters.doubles("drug_uptake_rate") = 0.01;
    parameters.doubles("drug_stress_threshold") = 0.03;
    parameters.doubles("efflux_induction_delay") = 30.0;
    parameters.doubles("efflux_strength") = 0.1;
    parameters.doubles("nrf2_decay_rate") = 0.01;
    parameters.doubles("drug_natural_decay_rate") = 0.0;
    parameters.doubles("hif1a_nrf2_priming_bonus") = 0.02;

    const double dt = 6.0;

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    Cell_Definition* pStroma = find_cell_definition("stromal_cell");
    assert(pTumor != NULL);
    assert(pStroma != NULL);

    Cell* tumor = create_cell(*pTumor);
    tumor->assign_position(std::vector<double>{100.0, 100.0, 0.0});

    const int hif1a_idx = tumor->custom_data.find_variable_index("hif1a_active");
    const int nrf2_idx = tumor->custom_data.find_variable_index("nrf2_active");
    const int abcb1_idx = tumor->custom_data.find_variable_index("abcb1_active");
    const int intra_idx = tumor->custom_data.find_variable_index("intracellular_drug");
    const int texp_idx = tumor->custom_data.find_variable_index("time_since_drug_exposure");
    assert(hif1a_idx >= 0);
    assert(nrf2_idx >= 0);
    assert(abcb1_idx >= 0);
    assert(intra_idx >= 0);
    assert(texp_idx >= 0);

    // Test A — No drug present.
    set_local_drug(tumor, 0.0);
    tumor->custom_data[intra_idx] = 0.0;
    tumor->custom_data[nrf2_idx] = 0.0;
    tumor->custom_data[abcb1_idx] = 0.0;
    tumor->custom_data[hif1a_idx] = 0.0;
    tumor->custom_data[texp_idx] = -1.0;
    module7_drug_response(tumor, tumor->phenotype, dt, ModulePhase::SENSING);
    assert(nearly_equal(tumor->custom_data[intra_idx], 0.0));
    assert(nearly_equal(tumor->custom_data[nrf2_idx], 0.0));
    assert(nearly_equal(tumor->custom_data[abcb1_idx], 0.0));
    assert(nearly_equal(tumor->custom_data[texp_idx], -1.0));
    std::cout << "PASS Test A" << std::endl;

    // Test B — Drug arrives, accumulates, and NRF2 ramps up under sustained exposure.
    set_local_drug(tumor, 0.5);
    tumor->custom_data[hif1a_idx] = 0.0;
    tumor->custom_data[intra_idx] = 0.0;
    tumor->custom_data[nrf2_idx] = 0.0;
    tumor->custom_data[abcb1_idx] = 0.0;
    tumor->custom_data[texp_idx] = -1.0;
    const double intra_before_b = tumor->custom_data[intra_idx];
    for (int i = 0; i < 10; ++i)
    {
        module7_drug_response(tumor, tumor->phenotype, dt, ModulePhase::SENSING);
    }
    assert(tumor->custom_data[intra_idx] > intra_before_b);
    assert(tumor->custom_data[texp_idx] > 0.0);
    assert(tumor->custom_data[nrf2_idx] > 0.4);
    std::cout << "PASS Test B" << std::endl;

    // Test C — ABCB1 induction is delayed and graded.
    Cell* tumor_c = create_cell(*pTumor);
    tumor_c->assign_position(std::vector<double>{140.0, 100.0, 0.0});
    const int c_hif1a = tumor_c->custom_data.find_variable_index("hif1a_active");
    const int c_nrf2 = tumor_c->custom_data.find_variable_index("nrf2_active");
    const int c_abcb1 = tumor_c->custom_data.find_variable_index("abcb1_active");
    const int c_intra = tumor_c->custom_data.find_variable_index("intracellular_drug");
    const int c_texp = tumor_c->custom_data.find_variable_index("time_since_drug_exposure");
    assert(c_hif1a >= 0);
    assert(c_nrf2 >= 0);
    assert(c_abcb1 >= 0);
    assert(c_intra >= 0);
    assert(c_texp >= 0);
    set_local_drug(tumor_c, 0.5);
    tumor_c->custom_data[c_hif1a] = 0.0;
    tumor_c->custom_data[c_nrf2] = 1.0;
    tumor_c->custom_data[c_abcb1] = 0.0;
    tumor_c->custom_data[c_intra] = 0.2;
    tumor_c->custom_data[c_texp] = 24.0; // below delay
    module7_drug_response(tumor_c, tumor_c->phenotype, dt, ModulePhase::SENSING);
    const double abcb1_pre_delay = tumor_c->custom_data[c_abcb1];
    assert(abcb1_pre_delay > 0.0);
    module7_drug_response(tumor_c, tumor_c->phenotype, dt, ModulePhase::SENSING);
    assert(tumor_c->custom_data[c_abcb1] > abcb1_pre_delay);
    std::cout << "PASS Test C" << std::endl;

    // Test D — Efflux reduces intracellular drug only while drug is present.
    set_local_drug(tumor_c, 0.5);
    tumor_c->custom_data[c_abcb1] = 1.0;
    tumor_c->custom_data[c_nrf2] = 1.0;
    tumor_c->custom_data[c_intra] = 0.8;
    tumor_c->custom_data[c_texp] = 36.0;

    Cell* tumor_no_efflux = create_cell(*pTumor);
    tumor_no_efflux->assign_position(std::vector<double>{180.0, 100.0, 0.0});
    const int ne_hif1a = tumor_no_efflux->custom_data.find_variable_index("hif1a_active");
    const int ne_nrf2 = tumor_no_efflux->custom_data.find_variable_index("nrf2_active");
    const int ne_abcb1 = tumor_no_efflux->custom_data.find_variable_index("abcb1_active");
    const int ne_intra = tumor_no_efflux->custom_data.find_variable_index("intracellular_drug");
    const int ne_texp = tumor_no_efflux->custom_data.find_variable_index("time_since_drug_exposure");
    assert(ne_hif1a >= 0);
    assert(ne_nrf2 >= 0);
    assert(ne_abcb1 >= 0);
    assert(ne_intra >= 0);
    assert(ne_texp >= 0);
    set_local_drug(tumor_no_efflux, 0.5);
    tumor_no_efflux->custom_data[ne_hif1a] = 0.0;
    tumor_no_efflux->custom_data[ne_nrf2] = 1.0;
    tumor_no_efflux->custom_data[ne_abcb1] = 0.0;
    tumor_no_efflux->custom_data[ne_intra] = 0.8;
    tumor_no_efflux->custom_data[ne_texp] = 36.0;

    for (int i = 0; i < 3; ++i)
    {
        module7_drug_response(tumor_c, tumor_c->phenotype, dt, ModulePhase::SENSING);
        module7_drug_response(tumor_no_efflux, tumor_no_efflux->phenotype, dt, ModulePhase::SENSING);
    }
    assert(tumor_c->custom_data[c_intra] < tumor_no_efflux->custom_data[ne_intra]);
    std::cout << "PASS Test D" << std::endl;

    // Test E — Stromal cell returns immediately.
    Cell* stroma = create_cell(*pStroma);
    stroma->assign_position(std::vector<double>{220.0, 100.0, 0.0});
    const double stroma_speed_before = stroma->phenotype.motility.migration_speed;
    set_local_drug(stroma, 0.8);
    module7_drug_response(stroma, stroma->phenotype, dt, ModulePhase::SENSING);
    assert(nearly_equal(stroma->phenotype.motility.migration_speed, stroma_speed_before));
    std::cout << "PASS Test E" << std::endl;

    // Test F — HIF1A priming lowers the effective activation requirement.
    Cell* tumor_f = create_cell(*pTumor);
    tumor_f->assign_position(std::vector<double>{260.0, 100.0, 0.0});
    const int f_hif1a = tumor_f->custom_data.find_variable_index("hif1a_active");
    const int f_nrf2 = tumor_f->custom_data.find_variable_index("nrf2_active");
    const int f_intra = tumor_f->custom_data.find_variable_index("intracellular_drug");
    const int f_texp = tumor_f->custom_data.find_variable_index("time_since_drug_exposure");
    assert(f_hif1a >= 0);
    assert(f_nrf2 >= 0);
    assert(f_intra >= 0);
    assert(f_texp >= 0);
    set_local_drug(tumor_f, 0.1);
    tumor_f->custom_data[f_intra] = 0.02;
    tumor_f->custom_data[f_nrf2] = 0.0;
    tumor_f->custom_data[f_hif1a] = 0.0;
    tumor_f->custom_data[f_texp] = -1.0;
    module7_drug_response(tumor_f, tumor_f->phenotype, dt, ModulePhase::SENSING);
    const double nrf2_without_hif = tumor_f->custom_data[f_nrf2];
    assert(nearly_equal(nrf2_without_hif, 0.0));
    tumor_f->custom_data[f_nrf2] = 0.0;
    tumor_f->custom_data[f_hif1a] = 1.0;
    module7_drug_response(tumor_f, tumor_f->phenotype, dt, ModulePhase::SENSING);
    assert(tumor_f->custom_data[f_nrf2] > nrf2_without_hif);
    assert(tumor_f->custom_data[f_nrf2] > 0.0);
    std::cout << "PASS Test F" << std::endl;

    std::cout << "PASS module7_drug_response_test" << std::endl;
    return 0;
}
