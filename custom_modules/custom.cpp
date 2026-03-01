#include "./custom.h"
#include "baseline_validation.h"
#include "tumor_calibration_knobs.h"

#include "../BioFVM/json.hpp"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace BioFVM;
using namespace PhysiCell;
using json = nlohmann::json;

// ---------------------------------------------------------------------------
// Global substrate indices (single source of truth for custom modules)
// ---------------------------------------------------------------------------
int oxygen_index = -1;
int tgfb_index   = -1;
int shh_index    = -1;
int ecm_index    = -1;
int drug_index   = -1;
ECMInterventionState intervention_state;

namespace
{

typedef void (*PhenotypeDispatchFn)(Cell*, Phenotype&, double);

std::unordered_map<int, PhenotypeDispatchFn> g_dispatch_by_type;
std::vector<double> g_ecm_ha_fraction_by_voxel;
std::vector<double> g_effective_o2_diffusion_by_voxel;
std::vector<double> g_effective_drug_diffusion_by_voxel;
std::vector<double> g_effective_tgfb_diffusion_by_voxel;
std::vector<double> g_effective_shh_diffusion_by_voxel;
std::vector<double> g_local_mechanical_stiffness_by_voxel;
bool g_module_trace_enabled = false;
std::vector<ModuleTraceEntry> g_module_trace_log;

// Store whichever diffusion-decay solver was configured before we wrapped it.
void (*g_base_diffusion_decay_solver)(Microenvironment&, double) = NULL;

TumorCalibrationProfiles g_knob_profiles;
bool g_knob_profiles_loaded = false;

inline double clamp_unit(double x)
{
    if (x < 0.0) return 0.0;
    if (x > 1.0) return 1.0;
    return x;
}

inline double clamp_nonnegative(double x)
{
    return (x < 0.0) ? 0.0 : x;
}

inline bool ends_with_json(const std::string& s)
{
    if (s.size() < 5) return false;
    return s.substr(s.size() - 5) == ".json";
}

void ensure_custom_scalar(Cell_Definition* pCD,
                          const std::string& name,
                          const std::string& units,
                          double default_value)
{
    if (pCD == NULL) return;
    if (pCD->custom_data.find_variable_index(name) < 0)
    {
        pCD->custom_data.add_variable(name, units, default_value);
    }
}

double read_custom_or_fallback(Cell* pCell, const std::string& name, double fallback)
{
    if (pCell == NULL) return fallback;
    const int idx = pCell->custom_data.find_variable_index(name);
    if (idx < 0) return fallback;
    return pCell->custom_data[idx];
}

void write_custom_if_present(Cell* pCell, const std::string& name, double value)
{
    if (pCell == NULL) return;
    const int idx = pCell->custom_data.find_variable_index(name);
    if (idx >= 0)
    {
        pCell->custom_data[idx] = value;
    }
}

void assign_tumor_identity_if_unset(Cell* pCell)
{
    if (pCell == NULL) return;

    const int cell_id_idx = pCell->custom_data.find_variable_index("cell_id");
    const int parent_id_idx = pCell->custom_data.find_variable_index("parent_id");
    if (cell_id_idx >= 0 && pCell->custom_data[cell_id_idx] < 0.0)
    {
        pCell->custom_data[cell_id_idx] = static_cast<double>(pCell->ID);
    }
    if (parent_id_idx >= 0 && pCell->custom_data[parent_id_idx] < 0.0)
    {
        pCell->custom_data[parent_id_idx] = -1.0;
    }
}

void assign_stromal_identity_if_unset(Cell* pCell)
{
    if (pCell == NULL) return;

    const int cell_id_idx = pCell->custom_data.find_variable_index("cell_id");
    const int parent_id_idx = pCell->custom_data.find_variable_index("parent_id");
    if (cell_id_idx >= 0 && pCell->custom_data[cell_id_idx] < 0.0)
    {
        pCell->custom_data[cell_id_idx] = static_cast<double>(pCell->ID);
    }
    if (parent_id_idx >= 0 && pCell->custom_data[parent_id_idx] < 0.0)
    {
        pCell->custom_data[parent_id_idx] = -1.0;
    }
}

void tumor_division_callback(Cell* pParent, Cell* pChild)
{
    if (pParent == NULL || pChild == NULL) return;

    const double parent_cell_id = [&]() -> double
    {
        const double stored = read_custom_or_fallback(pParent, "cell_id", -1.0);
        if (stored >= 0.0) return stored;
        return static_cast<double>(pParent->ID);
    }();

    write_custom_if_present(pParent, "cell_id", parent_cell_id);

    // Rule 4: daughter receives new unique ID and records parent lineage.
    write_custom_if_present(pChild, "cell_id", static_cast<double>(pChild->ID));
    write_custom_if_present(pChild, "parent_id", parent_cell_id);

    // Rule 5: reset all dynamic tumor state on daughter.
    write_custom_if_present(pChild, "hif1a_active", 0.0);
    write_custom_if_present(pChild, "zeb1_active", 0.0);
    write_custom_if_present(pChild, "cdh1_expressed", 1.0);
    write_custom_if_present(pChild, "nrf2_active", 0.0);
    write_custom_if_present(pChild, "abcb1_active", 0.0);
    write_custom_if_present(pChild, "mmp2_active", 0.0);
    write_custom_if_present(pChild, "intracellular_drug", 0.0);
    write_custom_if_present(pChild, "time_since_drug_exposure", -1.0);
    write_custom_if_present(pChild, "emt_extent", 0.0);
    write_custom_if_present(pChild, "time_alive", 0.0);
}

void stromal_division_callback(Cell* pParent, Cell* pChild)
{
    if (pParent == NULL || pChild == NULL) return;

    const double parent_cell_id = [&]() -> double
    {
        const double stored = read_custom_or_fallback(pParent, "cell_id", -1.0);
        if (stored >= 0.0) return stored;
        return static_cast<double>(pParent->ID);
    }();

    write_custom_if_present(pParent, "cell_id", parent_cell_id);
    write_custom_if_present(pChild, "cell_id", static_cast<double>(pChild->ID));
    write_custom_if_present(pChild, "parent_id", parent_cell_id);
}

void rebuild_stromal_custom_schema(Cell_Definition* pStroma, const Cell_Definition* pTumor)
{
    if (pStroma == NULL) return;

    Custom_Cell_Data schema;
    auto add_scalar = [&](const std::string& name, const std::string& units, double value)
    {
        schema.add_variable(name, units, value);
    };

    // Actor 2 identity.
    add_scalar("cell_id", "none", -1.0);
    add_scalar("parent_id", "none", -1.0);

    // Actor 2 stromal genotype subset (categorical codes, 0=WT).
    add_scalar("gene_state_SMAD4", "none", 0.0);
    add_scalar("gene_state_ACTA2", "none", 0.0);
    add_scalar("gene_state_GLI1", "none", 0.0);
    add_scalar("gene_state_HAS2", "none", 0.0);
    add_scalar("gene_state_COL1A1", "none", 0.0);
    add_scalar("gene_state_TGFB1", "none", 0.0);
    add_scalar("gene_state_SHH", "none", 0.0);

    // Activation / state variables.
    add_scalar("acta2_active", "dimensionless", 0.0);
    add_scalar("gli1_active", "dimensionless", 0.0);
    add_scalar("activation_mode", "dimensionless", 0.0);
    add_scalar("hif1a_active", "dimensionless", 0.0);
    add_scalar("ecm_production_rate", "dimensionless", 0.0);
    add_scalar("tgfb_secretion_active", "dimensionless", 0.0);
    add_scalar("time_alive", "min", 0.0);
    add_scalar("mechanical_pressure", "dimensionless", 0.0);

    // Keep vector lengths aligned with tumor custom data to avoid per-cell
    // serialization/index mismatches in output routines.
    const size_t target_size =
        (pTumor != NULL) ? pTumor->custom_data.variables.size() : schema.variables.size();
    int pad_counter = 0;
    while (schema.variables.size() < target_size)
    {
        add_scalar("_stromal_reserved_" + std::to_string(pad_counter), "none", 0.0);
        ++pad_counter;
    }

    pStroma->custom_data = schema;
}

std::vector<std::string> get_process_arguments()
{
    std::vector<std::string> args;

    // Linux-friendly way to recover argv[] from anywhere in the codebase.
    std::ifstream cmdline("/proc/self/cmdline", std::ios::binary);
    if (!cmdline.good())
    {
        return args;
    }

    std::string token;
    while (std::getline(cmdline, token, '\0'))
    {
        if (!token.empty())
        {
            args.push_back(token);
        }
    }

    return args;
}

std::string resolve_intervention_json_path()
{
    // First priority: explicit command-line options.
    const std::vector<std::string> args = get_process_arguments();
    for (size_t i = 1; i < args.size(); ++i)
    {
        const std::string& arg = args[i];

        if (arg == "--interventions" || arg == "--intervention-json")
        {
            if (i + 1 < args.size())
            {
                return args[i + 1];
            }
        }

        const std::string long_a = "--interventions=";
        const std::string long_b = "--intervention-json=";
        if (arg.compare(0, long_a.size(), long_a) == 0)
        {
            return arg.substr(long_a.size());
        }
        if (arg.compare(0, long_b.size(), long_b) == 0)
        {
            return arg.substr(long_b.size());
        }
    }

    // Second priority: positional arg after the XML config path.
    // Typical run is: ./program config/PhysiCell_settings.xml interventions.json
    if (args.size() > 2 && ends_with_json(args[2]))
    {
        return args[2];
    }

    // Third priority: XML user_parameters string fields (if present).
    if (parameters.strings.find_index("intervention_json") >= 0)
    {
        return parameters.strings("intervention_json");
    }
    if (parameters.strings.find_index("interventions_json") >= 0)
    {
        return parameters.strings("interventions_json");
    }

    return "";
}

void initialize_boolean_network_for_cell(Cell* pCell, CellType cell_type)
{
    if (pCell == NULL) return;

    BooleanNetwork* bn = get_boolean_network(pCell, cell_type);
    if (bn == NULL) return;

    // Explicitly initialize according to cell type, then write to custom_data.
    GeneParams params;
    bn->initialize(cell_type, params);
    bn->sync_to_cell(pCell);
}

void apply_xml_knob_overrides(TumorCalibrationKnobs& knobs)
{
    auto maybe_set_double = [&](const char* name, double& field)
    {
        try
        {
            field = parameters.doubles(name);
        }
        catch (...) {}
    };
    maybe_set_double("tgfb_secretion_rate", knobs.tgfb_secretion_rate);
    maybe_set_double("shh_secretion_rate", knobs.shh_secretion_rate);
    maybe_set_double("tgfb_brake_sensitivity", knobs.tgfb_brake_sensitivity);
    maybe_set_double("proliferation_rate", knobs.proliferation_rate);
    maybe_set_double("checkpoint_integrity", knobs.checkpoint_integrity);
    maybe_set_double("apoptosis_resistance", knobs.apoptosis_resistance);
    maybe_set_double("emt_induction_threshold", knobs.emt_induction_threshold);
    maybe_set_double("hypoxia_response_threshold", knobs.hypoxia_response_threshold);
    maybe_set_double("efflux_induction_delay", knobs.efflux_induction_delay);
    maybe_set_double("efflux_strength", knobs.efflux_strength);
    maybe_set_double("mechanical_compaction_strength", knobs.mechanical_compaction_strength);

    if (parameters.strings.find_index("emt_phenotype_extent") >= 0)
    {
        knobs.emt_phenotype_extent = knob_level_from_string(parameters.strings("emt_phenotype_extent"));
    }
    if (parameters.strings.find_index("hypoxia_phenotype_shift") >= 0)
    {
        knobs.hypoxia_phenotype_shift = knob_level_from_string(parameters.strings("hypoxia_phenotype_shift"));
    }
}

void load_knob_profiles_if_needed()
{
    if (g_knob_profiles_loaded) return;

    std::string profiles_path = "config/tumor_calibration_knobs.json";
    if (parameters.strings.find_index("tumor_calibration_knobs_json") >= 0)
    {
        profiles_path = parameters.strings("tumor_calibration_knobs_json");
    }
    g_knob_profiles = load_tumor_calibration_profiles(profiles_path);
    g_knob_profiles_loaded = true;
}

void load_interventions_at_setup()
{
    const std::string intervention_path = resolve_intervention_json_path();
    std::vector<KnobIntervention> knob_interventions;
    std::string profile_override;

    if (!intervention_path.empty())
    {
        std::ifstream ifs(intervention_path);
        if (!ifs.is_open())
        {
            std::cerr << "[setup_tissue] FATAL: cannot open intervention JSON: "
                      << intervention_path << "\n";
            exit(EXIT_FAILURE);
        }

        json root;
        ifs >> root;
        if (!root.is_object())
        {
            std::cerr << "[setup_tissue] FATAL: intervention JSON root must be an object.\n";
            exit(EXIT_FAILURE);
        }

        if (root.contains("calibration_profile"))
        {
            if (!root["calibration_profile"].is_string())
            {
                std::cerr << "[setup_tissue] FATAL: calibration_profile must be a string.\n";
                exit(EXIT_FAILURE);
            }
            profile_override = root["calibration_profile"].get<std::string>();
        }

        // Canonical schema.
        if (root.contains("knob_interventions"))
        {
            if (!root["knob_interventions"].is_array())
            {
                std::cerr << "[setup_tissue] FATAL: knob_interventions must be an array.\n";
                exit(EXIT_FAILURE);
            }
            for (const auto& entry : root["knob_interventions"])
            {
                if (!entry.is_object() ||
                    !entry.contains("knob") || !entry["knob"].is_string() ||
                    !entry.contains("effect") || !entry["effect"].is_string() ||
                    !entry.contains("strength") || !entry["strength"].is_number())
                {
                    std::cerr << "[setup_tissue] FATAL: invalid knob_interventions entry.\n";
                    exit(EXIT_FAILURE);
                }

                KnobIntervention iv;
                iv.knob = entry["knob"].get<std::string>();
                iv.effect = entry["effect"].get<std::string>();
                iv.strength = entry["strength"].get<double>();
                if (entry.contains("name") && entry["name"].is_string())
                {
                    iv.name = entry["name"].get<std::string>();
                }
                knob_interventions.push_back(iv);
            }
        }

        // One-release legacy bridge: interventions(gene) -> knob_interventions.
        if (root.contains("interventions"))
        {
            if (!root["interventions"].is_array())
            {
                std::cerr << "[setup_tissue] FATAL: legacy field 'interventions' must be an array.\n";
                exit(EXIT_FAILURE);
            }
            for (const auto& entry : root["interventions"])
            {
                if (!entry.is_object() ||
                    !entry.contains("gene") || !entry["gene"].is_string() ||
                    !entry.contains("effect") || !entry["effect"].is_string() ||
                    !entry.contains("strength") || !entry["strength"].is_number())
                {
                    continue;
                }

                const std::string gene = entry["gene"].get<std::string>();
                std::string mapped_knob;
                if (gene == "TGFB1") mapped_knob = "tgfb_secretion_rate";
                else if (gene == "SHH") mapped_knob = "shh_secretion_rate";
                else if (gene == "ABCB1") mapped_knob = "efflux_strength";
                else
                {
                    std::cerr << "[setup_tissue] FATAL: partition violation in legacy payload: gene '"
                              << gene
                              << "' is not bridgeable. Allowed legacy genes: TGFB1, SHH, ABCB1.\n";
                    exit(EXIT_FAILURE);
                }

                KnobIntervention iv;
                iv.knob = mapped_knob;
                iv.effect = entry["effect"].get<std::string>();
                iv.strength = entry["strength"].get<double>();
                if (entry.contains("name") && entry["name"].is_string())
                {
                    iv.name = entry["name"].get<std::string>();
                }
                knob_interventions.push_back(iv);
            }
        }
    }

    try
    {
        validate_knob_partition();
        load_knob_profiles_if_needed();

        std::string profile_name = "AsPC-1";
        if (parameters.strings.find_index("calibration_profile") >= 0)
        {
            profile_name = parameters.strings("calibration_profile");
        }
        if (!profile_override.empty())
        {
            profile_name = profile_override;
        }

        TumorCalibrationKnobs knobs = select_profile(g_knob_profiles, profile_name);
        apply_xml_knob_overrides(knobs);

        apply_targetable_knob_interventions(knobs, knob_interventions);
        set_active_calibration_knobs(knobs);

        // Load ECM component degradation parameters from XML user_parameters.
        // These are physical intervention parameters (Anchor 8 / validation only) that
        // bypass the knob partition — the EA never sets them, only the validation suite does.
        {
            TumorCalibrationKnobs ecm_knobs = get_active_calibration_knobs();
            if (parameters.doubles.find_index("ha_degrade_strength") >= 0)
                ecm_knobs.ha_degrade_strength  = parameters.doubles("ha_degrade_strength");
            if (parameters.doubles.find_index("col_degrade_strength") >= 0)
                ecm_knobs.col_degrade_strength = parameters.doubles("col_degrade_strength");
            set_active_calibration_knobs(ecm_knobs);
        }

        // Runtime interventions are knob-driven; keep gene-level interventions empty.
        set_current_interventions(std::vector<Intervention>());

        if (intervention_path.empty())
        {
            std::cerr << "[setup_tissue] No intervention JSON provided; using active calibration defaults.\n";
        }
        else
        {
            std::cerr << "[setup_tissue] Loaded knob interventions from: " << intervention_path
                      << " (count=" << knob_interventions.size() << ")\n";
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "[setup_tissue] FATAL: failed to initialize calibration knobs ("
                  << e.what() << ")\n";
        exit(EXIT_FAILURE);
    }
}

void apply_ecm_dependent_modifiers(Microenvironment& M, double dt)
{
    if (dt <= 0.0) return;

    static bool warned_missing_substrates = false;
    if (oxygen_index < 0 || drug_index < 0 || tgfb_index < 0 || shh_index < 0 || ecm_index < 0)
    {
        if (!warned_missing_substrates)
        {
            std::cerr << "[ecm_dependent_diffusion] WARNING: missing substrate indices "
                      << "(oxygen/tgfb/shh/drug/ecm_density). Skipping barrier modifier.\n";
            warned_missing_substrates = true;
        }
        return;
    }

    auto base_diffusion_or_zero = [&M](int substrate_index) -> double
    {
        if (substrate_index < 0) return 0.0;
        if (substrate_index >= static_cast<int>(M.diffusion_coefficients.size())) return 0.0;
        return std::max(0.0, M.diffusion_coefficients[substrate_index]);
    };

    const unsigned int n_voxels = M.number_of_voxels();
    for (unsigned int n = 0; n < n_voxels; ++n)
    {
        std::vector<double>& rho = M.density_vector(static_cast<int>(n));
        if (ecm_index >= static_cast<int>(rho.size()) ||
            drug_index >= static_cast<int>(rho.size()) ||
            oxygen_index >= static_cast<int>(rho.size()) ||
            tgfb_index >= static_cast<int>(rho.size()) ||
            shh_index >= static_cast<int>(rho.size()))
        {
            continue;
        }

        const double base_o2 = base_diffusion_or_zero(oxygen_index);
        const double base_tgfb = base_diffusion_or_zero(tgfb_index);
        const double base_shh = base_diffusion_or_zero(shh_index);
        const double base_drug = base_diffusion_or_zero(drug_index);

        if (base_o2 > 0.0 && n < g_effective_o2_diffusion_by_voxel.size())
        {
            const double ratio = clamp_unit(g_effective_o2_diffusion_by_voxel[n] / base_o2);
            rho[oxygen_index] = clamp_nonnegative(rho[oxygen_index] * ratio);
        }
        if (base_tgfb > 0.0 && n < g_effective_tgfb_diffusion_by_voxel.size())
        {
            const double ratio = clamp_unit(g_effective_tgfb_diffusion_by_voxel[n] / base_tgfb);
            rho[tgfb_index] = clamp_nonnegative(rho[tgfb_index] * ratio);
        }
        if (base_shh > 0.0 && n < g_effective_shh_diffusion_by_voxel.size())
        {
            const double ratio = clamp_unit(g_effective_shh_diffusion_by_voxel[n] / base_shh);
            rho[shh_index] = clamp_nonnegative(rho[shh_index] * ratio);
        }
        if (base_drug > 0.0 && n < g_effective_drug_diffusion_by_voxel.size())
        {
            const double ratio = clamp_unit(g_effective_drug_diffusion_by_voxel[n] / base_drug);
            rho[drug_index] = clamp_nonnegative(rho[drug_index] * ratio);
        }
    }
}

} // namespace

namespace
{

const char* phase_to_string(ModulePhase phase)
{
    switch (phase)
    {
        case ModulePhase::SENSING: return "SENSING";
        case ModulePhase::DECISION: return "DECISION";
        case ModulePhase::WRITE: return "WRITE";
        default: return "UNKNOWN";
    }
}

void assert_expected_phase(const char* module_name, ModulePhase actual, ModulePhase expected)
{
    if (actual != expected)
    {
        std::cerr << "[phase] FATAL: " << module_name
                  << " expected phase " << phase_to_string(expected)
                  << " but received " << phase_to_string(actual) << ".\n";
        assert(false && "Module invoked with wrong phase");
    }
}

void assert_microenvironment_write_allowed(ModulePhase phase, const char* module_name)
{
    if (phase != ModulePhase::WRITE)
    {
        std::cerr << "[phase] FATAL: " << module_name
                  << " attempted a microenvironment write during "
                  << phase_to_string(phase) << " phase.\n";
        assert(false && "Microenvironment writes are only allowed in WRITE phase");
    }
}

void log_module_trace(const char* module_name,
                      ModulePhase phase,
                      std::initializer_list<const char*> field_reads,
                      std::initializer_list<const char*> field_writes)
{
    if (!g_module_trace_enabled) return;

    ModuleTraceEntry entry;
    entry.module_name = module_name;
    entry.phase = phase;

    for (const char* field : field_reads)
    {
        if (field != NULL && field[0] != '\0')
        {
            entry.field_reads.push_back(field);
        }
    }

    for (const char* field : field_writes)
    {
        if (field != NULL && field[0] != '\0')
        {
            entry.field_writes.push_back(field);
        }
    }

    #pragma omp critical(module_trace_log_guard)
    {
        g_module_trace_log.push_back(entry);
    }
}

double read_density_value(const std::vector<double>& densities, int index)
{
    if (index < 0) return 0.0;
    if (index >= static_cast<int>(densities.size())) return 0.0;
    return densities[index];
}

void set_custom_data_if_present(Cell* pCell, const std::string& name, double value)
{
    if (pCell == NULL) return;
    const int idx = pCell->custom_data.find_variable_index(name);
    if (idx >= 0)
    {
        pCell->custom_data[idx] = value;
    }
}

double map_hypoxia_threshold_from_knob_6a(double knob_6a)
{
    return 0.1 * clamp_unit(knob_6a);
}

double read_custom_data_value_or_default(Cell* pCell,
                                         const std::string& name,
                                         double default_value)
{
    if (pCell == NULL) return default_value;
    const int idx = pCell->custom_data.find_variable_index(name);
    if (idx < 0) return default_value;
    return pCell->custom_data[idx];
}

double read_xml_double_or_default(const std::string& name, double default_value)
{
    if (parameters.doubles.find_index(name) >= 0)
    {
        return parameters.doubles(name);
    }
    return default_value;
}

void ensure_ecm_ha_fraction_storage(double default_fraction)
{
    const unsigned int n_voxels = microenvironment.number_of_voxels();
    const double clamped_default = clamp_unit(default_fraction);
    if (g_ecm_ha_fraction_by_voxel.size() != n_voxels)
    {
        g_ecm_ha_fraction_by_voxel.assign(n_voxels, clamped_default);
    }
}

int get_cell_voxel_index(Cell* pCell)
{
    if (pCell == NULL) return -1;
    std::vector<double> position = pCell->position;
    return microenvironment.nearest_voxel_index(position);
}

} // namespace

void set_module_trace_enabled(bool enabled)
{
    g_module_trace_enabled = enabled;
}

void clear_module_trace_log(void)
{
    #pragma omp critical(module_trace_log_guard)
    {
        g_module_trace_log.clear();
    }
}

std::vector<ModuleTraceEntry> get_module_trace_log(void)
{
    std::vector<ModuleTraceEntry> copy;
    #pragma omp critical(module_trace_log_guard)
    {
        copy = g_module_trace_log;
    }
    return copy;
}

void reset_ecm_ha_fraction_field(double default_fraction)
{
    const unsigned int n_voxels = microenvironment.number_of_voxels();
    const double clamped_default = std::max(0.0, std::min(1.0, default_fraction));
    g_ecm_ha_fraction_by_voxel.assign(n_voxels, clamped_default);
}

double get_ecm_ha_fraction(int voxel_index)
{
    if (voxel_index < 0) return 0.0;
    if (voxel_index >= static_cast<int>(g_ecm_ha_fraction_by_voxel.size())) return 0.0;
    return g_ecm_ha_fraction_by_voxel[voxel_index];
}

void set_ecm_ha_fraction(int voxel_index, double value)
{
    if (voxel_index < 0) return;
    if (voxel_index >= static_cast<int>(g_ecm_ha_fraction_by_voxel.size())) return;
    g_ecm_ha_fraction_by_voxel[voxel_index] = std::max(0.0, std::min(1.0, value));
}

double get_ecm_collagen_fraction(int voxel_index)
{
    return 1.0 - get_ecm_ha_fraction(voxel_index);
}

void update_ecm_effective_diffusion_coefficients(Microenvironment& M)
{
    const unsigned int n_voxels = M.number_of_voxels();

    auto base_diffusion_or_zero = [&M](int substrate_index) -> double
    {
        if (substrate_index < 0) return 0.0;
        if (substrate_index >= static_cast<int>(M.diffusion_coefficients.size())) return 0.0;
        return std::max(0.0, M.diffusion_coefficients[substrate_index]);
    };

    const double base_o2 = base_diffusion_or_zero(oxygen_index);
    const double base_tgfb = base_diffusion_or_zero(tgfb_index);
    const double base_shh = base_diffusion_or_zero(shh_index);
    const double base_drug = base_diffusion_or_zero(drug_index);

    g_effective_o2_diffusion_by_voxel.assign(n_voxels, base_o2);
    g_effective_tgfb_diffusion_by_voxel.assign(n_voxels, base_tgfb);
    g_effective_shh_diffusion_by_voxel.assign(n_voxels, base_shh);
    g_effective_drug_diffusion_by_voxel.assign(n_voxels, base_drug);
    g_local_mechanical_stiffness_by_voxel.assign(n_voxels, 0.0);

    if (ecm_index < 0)
    {
        return;
    }

    for (unsigned int n = 0; n < n_voxels; ++n)
    {
        const std::vector<double>& rho = M.density_vector(static_cast<int>(n));
        if (ecm_index >= static_cast<int>(rho.size()))
        {
            continue;
        }

        const double ecm = clamp_unit(rho[ecm_index]);
        const double ha_frac =
            (n < g_ecm_ha_fraction_by_voxel.size()) ? clamp_unit(g_ecm_ha_fraction_by_voxel[n]) : 0.5;
        const double col_frac = 1.0 - ha_frac;

        // Collagen-dominant stiffness term used by mechanical modules.
        g_local_mechanical_stiffness_by_voxel[n] = ecm * col_frac;

        if (base_o2 > 0.0)
        {
            const double weight_o2 = 0.6 * ha_frac + 0.3 * col_frac;
            double effective = base_o2 * (1.0 - ecm * weight_o2);
            effective = std::max(effective, base_o2 * 0.05);
            g_effective_o2_diffusion_by_voxel[n] = effective;
        }

        if (base_drug > 0.0)
        {
            const double weight_drug = 0.8 * ha_frac + 0.3 * col_frac;
            double effective = base_drug * (1.0 - ecm * weight_drug);
            effective = std::max(effective, base_drug * 0.02);
            g_effective_drug_diffusion_by_voxel[n] = effective;
        }

        if (base_tgfb > 0.0)
        {
            const double weight_tgfb = 0.7 * ha_frac + 0.4 * col_frac;
            double effective = base_tgfb * (1.0 - ecm * weight_tgfb);
            effective = std::max(effective, base_tgfb * 0.05);
            g_effective_tgfb_diffusion_by_voxel[n] = effective;
        }

        if (base_shh > 0.0)
        {
            const double weight_shh = 0.7 * ha_frac + 0.4 * col_frac;
            double effective = base_shh * (1.0 - ecm * weight_shh);
            effective = std::max(effective, base_shh * 0.05);
            g_effective_shh_diffusion_by_voxel[n] = effective;
        }
    }
}

double get_effective_diffusion_coefficient(int substrate_index, int voxel_index)
{
    if (voxel_index < 0) return 0.0;

    if (substrate_index == oxygen_index &&
        voxel_index < static_cast<int>(g_effective_o2_diffusion_by_voxel.size()))
    {
        return g_effective_o2_diffusion_by_voxel[voxel_index];
    }
    if (substrate_index == tgfb_index &&
        voxel_index < static_cast<int>(g_effective_tgfb_diffusion_by_voxel.size()))
    {
        return g_effective_tgfb_diffusion_by_voxel[voxel_index];
    }
    if (substrate_index == shh_index &&
        voxel_index < static_cast<int>(g_effective_shh_diffusion_by_voxel.size()))
    {
        return g_effective_shh_diffusion_by_voxel[voxel_index];
    }
    if (substrate_index == drug_index &&
        voxel_index < static_cast<int>(g_effective_drug_diffusion_by_voxel.size()))
    {
        return g_effective_drug_diffusion_by_voxel[voxel_index];
    }

    if (substrate_index >= 0 &&
        substrate_index < static_cast<int>(microenvironment.diffusion_coefficients.size()))
    {
        return std::max(0.0, microenvironment.diffusion_coefficients[substrate_index]);
    }

    return 0.0;
}

double get_local_mechanical_stiffness(int voxel_index)
{
    if (voxel_index < 0) return 0.0;
    if (voxel_index >= static_cast<int>(g_local_mechanical_stiffness_by_voxel.size())) return 0.0;
    return std::max(0.0, g_local_mechanical_stiffness_by_voxel[voxel_index]);
}

// SENSING PHASE (reads environment, writes nothing)
void module1_oxygen_sensing(Cell* pCell, Phenotype& phenotype, double dt, ModulePhase phase)
{
    (void)phenotype;
    (void)dt;
    assert_expected_phase("module1_oxygen_sensing", phase, ModulePhase::SENSING);
    log_module_trace("module1_oxygen_sensing", phase, {"oxygen"}, {});

    if (pCell == NULL) return;

    const std::vector<double>& densities = pCell->nearest_density_vector();
    const double local_o2 = read_density_value(densities, oxygen_index);

    double knob_6a = 0.0;
    if (parameters.doubles.find_index("hypoxia_response_threshold") >= 0)
    {
        knob_6a = parameters.doubles("hypoxia_response_threshold");
    }
    const double threshold_mapped_from_knob_6a = map_hypoxia_threshold_from_knob_6a(knob_6a);

    if (local_o2 < threshold_mapped_from_knob_6a)
    {
        set_custom_data_if_present(pCell, "hif1a_active", 1.0);
    }
    else
    {
        set_custom_data_if_present(pCell, "hif1a_active", 0.0);
    }
}

void module3_stromal_activation(Cell* pCell, Phenotype& phenotype, double dt, ModulePhase phase)
{
    (void)dt;
    assert_expected_phase("module3_stromal_activation", phase, ModulePhase::SENSING);
    log_module_trace("module3_stromal_activation", phase, {"tgfb", "shh"}, {});

    if (pCell == NULL) return;

    Cell_Definition* pStromaDef = find_cell_definition("stromal_cell");
    const bool is_stromal =
        (pStromaDef != NULL) && (pCell->type == pStromaDef->type);
    if (!is_stromal)
    {
        return;
    }

    const std::vector<double>& densities = pCell->nearest_density_vector();
    const double local_tgfb = read_density_value(densities, tgfb_index);
    const double local_shh = read_density_value(densities, shh_index);

    const double activation_threshold =
        read_xml_double_or_default("tgfb_activation_threshold", 0.0);
    const double shh_threshold =
        read_xml_double_or_default("shh_activation_threshold", 0.0);
    const double ecm_production_rate_base =
        read_xml_double_or_default("ecm_production_rate_base", 0.0);
    const double ecm_production_rate_boosted =
        read_xml_double_or_default("ecm_production_rate_boosted", 0.0);
    const double psc_migration_speed =
        read_xml_double_or_default("psc_migration_speed", 0.05);
    const double caf_migration_speed =
        read_xml_double_or_default("caf_migration_speed", 0.1);

    const double previous_acta2 =
        read_custom_data_value_or_default(pCell, "acta2_active", 0.0);
    double updated_gli1 =
        read_custom_data_value_or_default(pCell, "gli1_active", 0.0);
    double updated_acta2 = previous_acta2;

    if (is_stromal && previous_acta2 == 0.0)
    {
        const double combined_signal = local_tgfb + local_shh;
        if (combined_signal > activation_threshold)
        {
            updated_acta2 = 1.0; // IRREVERSIBLE

            if (local_shh > 0.0)
            {
                updated_gli1 = 1.0;
            }
        }
    }

    if (updated_acta2 == 1.0)
    {
        updated_acta2 = 1.0; // acta2 stays ON permanently

        if (local_shh > shh_threshold)
        {
            updated_gli1 = 1.0;
        }
        else
        {
            updated_gli1 = 0.0; // GLI1 is reversible
        }
    }

    if (previous_acta2 == 1.0 && updated_acta2 == 0.0)
    {
        std::cerr << "[module3_stromal_activation] FATAL: acta2_active attempted 1->0 transition.\n";
        assert(false && "acta2_active is irreversible and cannot transition from 1.0 to 0.0");
    }

    set_custom_data_if_present(pCell, "acta2_active", updated_acta2);
    set_custom_data_if_present(pCell, "gli1_active", updated_gli1);

    if (updated_acta2 == 1.0)
    {
        set_custom_data_if_present(pCell, "activation_mode", 1.0);
        set_custom_data_if_present(pCell, "tgfb_secretion_active", 1.0);
        if (updated_gli1 == 1.0)
        {
            set_custom_data_if_present(pCell, "ecm_production_rate", ecm_production_rate_boosted);
        }
        else
        {
            set_custom_data_if_present(pCell, "ecm_production_rate", ecm_production_rate_base);
        }
        phenotype.motility.migration_speed = caf_migration_speed;
    }
    else
    {
        set_custom_data_if_present(pCell, "activation_mode", 0.0);
        set_custom_data_if_present(pCell, "tgfb_secretion_active", 0.0);
        set_custom_data_if_present(pCell, "ecm_production_rate", 0.0);
        phenotype.motility.migration_speed = psc_migration_speed;
    }
}

void module7_drug_response(Cell* pCell, Phenotype& phenotype, double dt, ModulePhase phase)
{
    assert_expected_phase("module7_drug_response", phase, ModulePhase::SENSING);
    log_module_trace("module7_drug_response", phase, {"drug"}, {});

    if (pCell == NULL) return;

    Cell_Definition* pTumorDef = find_cell_definition("tumor_cell");
    const bool is_tumor =
        (pTumorDef != NULL) && (pCell->type == pTumorDef->type);
    if (!is_tumor)
    {
        return;
    }

    if (drug_index < 0) return;

    const std::vector<double>& densities = pCell->nearest_density_vector();
    const double local_drug = read_density_value(densities, drug_index);

    const double hif1a_active =
        read_custom_data_value_or_default(pCell, "hif1a_active", 0.0);
    double nrf2_active =
        read_custom_data_value_or_default(pCell, "nrf2_active", 0.0);
    double abcb1_active =
        read_custom_data_value_or_default(pCell, "abcb1_active", 0.0);
    double intracellular_drug =
        read_custom_data_value_or_default(pCell, "intracellular_drug", 0.0);
    double time_since_drug_exposure =
        read_custom_data_value_or_default(pCell, "time_since_drug_exposure", -1.0);

    const double drug_uptake_rate =
        read_xml_double_or_default("drug_uptake_rate", 0.0);
    const double drug_stress_threshold =
        read_xml_double_or_default("drug_stress_threshold", 0.0);
    const double efflux_induction_delay =
        read_xml_double_or_default("efflux_induction_delay", 0.0);
    const double efflux_strength =
        read_xml_double_or_default("efflux_strength", 0.0);
    const double nrf2_decay_rate =
        read_xml_double_or_default("nrf2_decay_rate", 0.0);
    const double drug_natural_decay_rate =
        read_xml_double_or_default("drug_natural_decay_rate", 0.0);
    const double hif1a_nrf2_priming_bonus =
        read_xml_double_or_default("hif1a_nrf2_priming_bonus", 0.0);

    // Step 1: Drug uptake.
    if (local_drug > 0.0)
    {
        const double uptake = drug_uptake_rate * local_drug * dt;
        intracellular_drug = intracellular_drug + uptake;
        intracellular_drug = std::min(1.0, intracellular_drug);

        if (time_since_drug_exposure < 0.0)
        {
            time_since_drug_exposure = 0.0;
        }
    }

    // Step 2: Intracellular drug decay.
    intracellular_drug = intracellular_drug * (1.0 - drug_natural_decay_rate * dt);
    intracellular_drug = std::max(0.0, intracellular_drug);

    // Step 3: Stress sensing (NRF2 activation).
    double stress_signal = intracellular_drug;
    if (hif1a_active == 1.0)
    {
        stress_signal = stress_signal + hif1a_nrf2_priming_bonus;
    }

    if (stress_signal > drug_stress_threshold)
    {
        nrf2_active = 1.0;
    }
    else
    {
        if (nrf2_active > 0.0)
        {
            nrf2_active = nrf2_active - nrf2_decay_rate * dt;
            nrf2_active = std::max(0.0, nrf2_active);
            if (nrf2_active < 0.5)
            {
                nrf2_active = 0.0;
            }
        }
    }

    // Step 4: Efflux activation (ABCB1).
    if (nrf2_active >= 1.0 && time_since_drug_exposure >= 0.0)
    {
        if (time_since_drug_exposure >= efflux_induction_delay)
        {
            abcb1_active = 1.0;
        }
    }

    // Step 5: Efflux pump action.
    if (abcb1_active == 1.0)
    {
        const double efflux_amount = efflux_strength * intracellular_drug * dt;
        intracellular_drug = intracellular_drug - efflux_amount;
        intracellular_drug = std::max(0.0, intracellular_drug);
    }

    // Step 6: Advance exposure clock.
    if (time_since_drug_exposure >= 0.0)
    {
        time_since_drug_exposure = time_since_drug_exposure + dt;
    }

    // Step 7: Drug uptake from environment (sink effect via PhysiCell uptake).
    if (drug_index >= 0 &&
        drug_index < static_cast<int>(phenotype.secretion.uptake_rates.size()))
    {
        if (local_drug > 0.0)
        {
            phenotype.secretion.uptake_rates[drug_index] = drug_uptake_rate;
        }
        else
        {
            phenotype.secretion.uptake_rates[drug_index] = 0.0;
        }
    }

    intracellular_drug = std::min(1.0, std::max(0.0, intracellular_drug));

    set_custom_data_if_present(pCell, "intracellular_drug", intracellular_drug);
    set_custom_data_if_present(pCell, "nrf2_active", nrf2_active);
    set_custom_data_if_present(pCell, "abcb1_active", abcb1_active);
    set_custom_data_if_present(pCell, "time_since_drug_exposure", time_since_drug_exposure);
}

// DECISION PHASE (reads module outputs + environment, writes nothing to fields)
void module5_emt_engine(Cell* pCell, Phenotype& phenotype, double dt, ModulePhase phase)
{
    (void)dt;
    assert_expected_phase("module5_emt_engine", phase, ModulePhase::DECISION);
    log_module_trace("module5_emt_engine", phase, {"tgfb"}, {});

    if (pCell == NULL) return;

    Cell_Definition* pTumorDef = find_cell_definition("tumor_cell");
    const bool is_tumor =
        (pTumorDef != NULL) && (pCell->type == pTumorDef->type);
    if (!is_tumor)
    {
        return;
    }

    const std::vector<double>& densities = pCell->nearest_density_vector();
    const double local_tgfb = read_density_value(densities, tgfb_index);
    const double hif1a_active =
        read_custom_data_value_or_default(pCell, "hif1a_active", 0.0);

    const double emt_induction_threshold =
        read_xml_double_or_default("emt_induction_threshold", 0.0);
    int emt_phenotype_extent = 1;
    if (parameters.ints.find_index("emt_phenotype_extent") >= 0)
    {
        emt_phenotype_extent = parameters.ints("emt_phenotype_extent");
    }
    const double hif1a_emt_boost =
        read_xml_double_or_default("hif1a_emt_boost", 0.0);
    const double motility_epithelial =
        read_xml_double_or_default("motility_epithelial", 0.0);
    const double motility_mesenchymal_low =
        read_xml_double_or_default("motility_mesenchymal_low", 0.0);
    const double motility_mesenchymal_med =
        read_xml_double_or_default("motility_mesenchymal_med", 0.0);
    const double motility_mesenchymal_high =
        read_xml_double_or_default("motility_mesenchymal_high", 0.0);
    const double adhesion_epithelial =
        read_xml_double_or_default("adhesion_epithelial", 0.0);
    const double adhesion_mesenchymal_low =
        read_xml_double_or_default("adhesion_mesenchymal_low", 0.0);
    const double adhesion_mesenchymal_med =
        read_xml_double_or_default("adhesion_mesenchymal_med", 0.0);
    const double adhesion_mesenchymal_high =
        read_xml_double_or_default("adhesion_mesenchymal_high", 0.0);

    double induction_signal = local_tgfb;
    if (hif1a_active == 1.0)
    {
        induction_signal = induction_signal + hif1a_emt_boost;
    }

    if (induction_signal > emt_induction_threshold)
    {
        set_custom_data_if_present(pCell, "zeb1_active", 1.0);
        set_custom_data_if_present(pCell, "cdh1_expressed", 0.0);

        if (emt_phenotype_extent == 1)
        {
            phenotype.motility.migration_speed = motility_mesenchymal_low;
            phenotype.mechanics.cell_cell_adhesion_strength = adhesion_mesenchymal_low;
            set_custom_data_if_present(pCell, "mmp2_active", 0.0);
            set_custom_data_if_present(pCell, "emt_extent", 1.0);
        }

        if (emt_phenotype_extent == 2)
        {
            phenotype.motility.migration_speed = motility_mesenchymal_med;
            phenotype.mechanics.cell_cell_adhesion_strength = adhesion_mesenchymal_med;
            set_custom_data_if_present(pCell, "mmp2_active", 1.0);
            set_custom_data_if_present(pCell, "emt_extent", 2.0);
        }

        if (emt_phenotype_extent == 3)
        {
            phenotype.motility.migration_speed = motility_mesenchymal_high;
            phenotype.mechanics.cell_cell_adhesion_strength = adhesion_mesenchymal_high;
            set_custom_data_if_present(pCell, "mmp2_active", 1.0);
            set_custom_data_if_present(pCell, "emt_extent", 3.0);
        }
    }
    else
    {
        set_custom_data_if_present(pCell, "zeb1_active", 0.0);
        set_custom_data_if_present(pCell, "cdh1_expressed", 1.0);
        set_custom_data_if_present(pCell, "mmp2_active", 0.0);
        set_custom_data_if_present(pCell, "emt_extent", 0.0);
        phenotype.motility.migration_speed = motility_epithelial;
        phenotype.mechanics.cell_cell_adhesion_strength = adhesion_epithelial;
    }
}

void module4_proliferation_death(Cell* pCell, Phenotype& phenotype, double dt, ModulePhase phase)
{
    (void)dt;
    assert_expected_phase("module4_proliferation_death", phase, ModulePhase::DECISION);
    log_module_trace("module4_proliferation_death", phase, {"tgfb"}, {});

    if (pCell == NULL) return;

    const std::vector<double>& densities = pCell->nearest_density_vector();
    const double local_tgfb = read_density_value(densities, tgfb_index);

    const double hif1a_active =
        read_custom_data_value_or_default(pCell, "hif1a_active", 0.0);
    const double zeb1_active =
        read_custom_data_value_or_default(pCell, "zeb1_active", 0.0);
    const double acta2_active =
        read_custom_data_value_or_default(pCell, "acta2_active", 0.0);
    const double gli1_active =
        read_custom_data_value_or_default(pCell, "gli1_active", 0.0);
    const double abcb1_active =
        read_custom_data_value_or_default(pCell, "abcb1_active", 0.0);
    const double intracellular_drug =
        read_custom_data_value_or_default(pCell, "intracellular_drug", 0.0);
    const double mechanical_pressure =
        read_custom_data_value_or_default(pCell, "mechanical_pressure", 0.0);
    const double smad4_state =
        read_custom_data_value_or_default(pCell, "SMAD4", 0.0);

    const int apoptosis_index =
        phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);
    if (apoptosis_index < 0 ||
        apoptosis_index >= static_cast<int>(phenotype.death.rates.size()))
    {
        return;
    }

    Cell_Definition* pTumorDef = find_cell_definition("tumor_cell");
    Cell_Definition* pStromaDef = find_cell_definition("stromal_cell");

    const bool is_tumor =
        (pTumorDef != NULL) && (pCell->type == pTumorDef->type);
    const bool is_stromal =
        (pStromaDef != NULL) && (pCell->type == pStromaDef->type);

    if (is_tumor)
    {
        const double base_proliferation_rate =
            read_xml_double_or_default("base_proliferation_rate", 0.0);
        const double apoptosis_resistance =
            read_xml_double_or_default("apoptosis_resistance", 0.0);
        const double go_grow_penalty =
            read_xml_double_or_default("go_grow_penalty", 0.0);
        const double hypoxia_proliferation_modifier =
            read_xml_double_or_default("hypoxia_proliferation_modifier", 0.0);
        const double hypoxia_death_resistance_bonus =
            read_xml_double_or_default("hypoxia_death_resistance_bonus", 0.0);
        const double contact_inhibition_threshold =
            read_xml_double_or_default("contact_inhibition_threshold", 0.0);
        const double drug_kill_coefficient =
            read_xml_double_or_default("drug_kill_coefficient", 0.0);
        const double efflux_drug_reduction =
            read_xml_double_or_default("efflux_drug_reduction", 0.0);
        const double tgfb_brake_sensitivity =
            read_xml_double_or_default("tgfb_brake_sensitivity", 0.0);

        const double gene_state_kras =
            read_custom_data_value_or_default(pCell, "gene_state_KRAS", 1.0);
        const double gene_state_myc =
            read_custom_data_value_or_default(pCell, "gene_state_MYC", 0.0);

        const bool kras_is_wt = (gene_state_kras < 0.5);
        const bool kras_is_gof = (gene_state_kras >= 0.5 && gene_state_kras < 1.5);
        const bool myc_is_wt = (gene_state_myc < 0.5);
        const bool myc_is_gof = (gene_state_myc >= 0.5 && gene_state_myc < 1.5);

        // Rule 1: KRAS/MYC genotype modulates growth drive.
        // Keep canonical PDAC baseline (KRAS=GOF, MYC=WT) at multiplier 1.0.
        double genotype_proliferation_multiplier = 1.0;
        if (kras_is_wt && myc_is_wt)
        {
            genotype_proliferation_multiplier = 0.5;
        }
        else if (kras_is_gof && myc_is_gof)
        {
            genotype_proliferation_multiplier = 1.5;
        }

        double proliferation = base_proliferation_rate;
        proliferation = proliferation * genotype_proliferation_multiplier;

        if (zeb1_active == 1.0)
        {
            proliferation = proliferation * (1.0 - go_grow_penalty);
        }

        if (hif1a_active == 1.0)
        {
            proliferation = proliferation * (1.0 - hypoxia_proliferation_modifier);
        }

        // E09: TGF-beta growth arrest is SMAD4-dependent and only active in SMAD4-WT.
        if (smad4_state > 0.5)
        {
            const double tgfb_brake = clamp_unit(std::max(0.0, local_tgfb) * tgfb_brake_sensitivity);
            proliferation = proliferation * (1.0 - tgfb_brake);
        }

        if (mechanical_pressure > contact_inhibition_threshold)
        {
            proliferation = 0.0;
        }

        phenotype.cycle.data.transition_rate(0, 0) = proliferation;

        double effective_drug = intracellular_drug;
        if (abcb1_active == 1.0)
        {
            effective_drug = effective_drug * (1.0 - efflux_drug_reduction);
        }

        const double base_death_rate = (1.0 - apoptosis_resistance) * 0.001;

        double death_resistance_modifier = 0.0;
        if (hif1a_active == 1.0)
        {
            death_resistance_modifier = hypoxia_death_resistance_bonus;
        }

        const double drug_kill_rate = effective_drug * drug_kill_coefficient;

        const double total_apoptosis_rate = std::max(
            0.0,
            base_death_rate * (1.0 - death_resistance_modifier) + drug_kill_rate);

        phenotype.death.rates[apoptosis_index] = total_apoptosis_rate;
    }

    if (is_stromal)
    {
        const double caf_proliferation_rate =
            read_xml_double_or_default("caf_proliferation_rate", 0.0);
        const double gli1_proliferation_boost =
            read_xml_double_or_default("gli1_proliferation_boost", 1.0);
        const double psc_proliferation_rate =
            read_xml_double_or_default("psc_proliferation_rate", 0.0);

        if (acta2_active == 1.0)
        {
            double caf_rate = caf_proliferation_rate;
            if (gli1_active == 1.0)
            {
                caf_rate = caf_rate * gli1_proliferation_boost;
            }
            phenotype.cycle.data.transition_rate(0, 0) = caf_rate;
            phenotype.death.rates[apoptosis_index] = 0.0;
        }
        else
        {
            phenotype.cycle.data.transition_rate(0, 0) = psc_proliferation_rate;
            phenotype.death.rates[apoptosis_index] = 0.0;
        }
    }
}

// WRITE PHASE (writes to microenvironment fields)
void module2_paracrine_secretion(Cell* pCell, Phenotype& phenotype, double dt, ModulePhase phase)
{
    (void)dt;
    assert_expected_phase("module2_paracrine_secretion", phase, ModulePhase::WRITE);
    assert_microenvironment_write_allowed(phase, "module2_paracrine_secretion");
    log_module_trace("module2_paracrine_secretion", phase, {}, {"tgfb", "shh"});

    if (pCell == NULL) return;

    Cell_Definition* pTumorDef = find_cell_definition("tumor_cell");
    Cell_Definition* pStromaDef = find_cell_definition("stromal_cell");

    const bool is_tumor =
        (pTumorDef != NULL) && (pCell->type == pTumorDef->type);
    const bool is_stromal =
        (pStromaDef != NULL) && (pCell->type == pStromaDef->type);

    const double hif1a_active =
        read_custom_data_value_or_default(pCell, "hif1a_active", 0.0);
    const double acta2_active =
        read_custom_data_value_or_default(pCell, "acta2_active", 0.0);

    const double tgfb_secretion_rate =
        read_xml_double_or_default("tgfb_secretion_rate", 0.0);
    const double shh_secretion_rate =
        read_xml_double_or_default("shh_secretion_rate", 0.0);
    const double hif1a_tgfb_amplification_factor =
        read_xml_double_or_default("hif1a_tgfb_amplification_factor", 1.0);
    const double caf_tgfb_secretion_rate =
        read_xml_double_or_default("caf_tgfb_secretion_rate", 0.0);

    if (is_tumor)
    {
        double tgfb_out = tgfb_secretion_rate;
        double shh_out  = shh_secretion_rate;
        if (hif1a_active == 1.0)
        {
            tgfb_out = tgfb_out * hif1a_tgfb_amplification_factor;
        }

        if (tgfb_index >= 0 &&
            tgfb_index < static_cast<int>(phenotype.secretion.secretion_rates.size()))
        {
            phenotype.secretion.secretion_rates[tgfb_index] = tgfb_out;
        }
        if (shh_index >= 0 &&
            shh_index < static_cast<int>(phenotype.secretion.secretion_rates.size()))
        {
            phenotype.secretion.secretion_rates[shh_index] = shh_out;
        }
    }

    if (is_stromal && acta2_active == 1.0)
    {
        if (tgfb_index >= 0 &&
            tgfb_index < static_cast<int>(phenotype.secretion.secretion_rates.size()))
        {
            phenotype.secretion.secretion_rates[tgfb_index] = caf_tgfb_secretion_rate;
        }
        if (shh_index >= 0 &&
            shh_index < static_cast<int>(phenotype.secretion.secretion_rates.size()))
        {
            phenotype.secretion.secretion_rates[shh_index] = 0.0;
        }
    }

    if (is_stromal && acta2_active == 0.0)
    {
        if (tgfb_index >= 0 &&
            tgfb_index < static_cast<int>(phenotype.secretion.secretion_rates.size()))
        {
            phenotype.secretion.secretion_rates[tgfb_index] = 0.0;
        }
        if (shh_index >= 0 &&
            shh_index < static_cast<int>(phenotype.secretion.secretion_rates.size()))
        {
            phenotype.secretion.secretion_rates[shh_index] = 0.0;
        }
    }
}

void module6_ecm_production(Cell* pCell, Phenotype& phenotype, double dt, ModulePhase phase)
{
    (void)phenotype;
    assert_expected_phase("module6_ecm_production", phase, ModulePhase::WRITE);
    assert_microenvironment_write_allowed(phase, "module6_ecm_production");
    log_module_trace("module6_ecm_production", phase, {"ecm_density", "ecm_ha_fraction"}, {"ecm_density", "ecm_ha_fraction"});

    if (pCell == NULL) return;
    if (ecm_index < 0) return;

    const double ecm_ha_fraction_default =
        read_xml_double_or_default("ecm_ha_fraction_default", 0.5);
    ensure_ecm_ha_fraction_storage(ecm_ha_fraction_default);

    const int voxel_index = get_cell_voxel_index(pCell);
    if (voxel_index < 0 ||
        voxel_index >= static_cast<int>(microenvironment.number_of_voxels()))
    {
        return;
    }

    std::vector<double>& rho = microenvironment.density_vector(voxel_index);
    if (ecm_index >= static_cast<int>(rho.size())) return;

    double ecm_density = clamp_unit(rho[ecm_index]);
    double ecm_ha_fraction = clamp_unit(g_ecm_ha_fraction_by_voxel[voxel_index]);

    const double acta2_active =
        read_custom_data_value_or_default(pCell, "acta2_active", 0.0);
    const double gli1_active =
        read_custom_data_value_or_default(pCell, "gli1_active", 0.0);
    const double mmp2_active =
        read_custom_data_value_or_default(pCell, "mmp2_active", 0.0);

    Cell_Definition* pTumorDef = find_cell_definition("tumor_cell");
    Cell_Definition* pStromaDef = find_cell_definition("stromal_cell");
    const bool is_tumor =
        (pTumorDef != NULL) && (pCell->type == pTumorDef->type);
    const bool is_stromal =
        (pStromaDef != NULL) && (pCell->type == pStromaDef->type);

    const double ecm_production_rate_base =
        read_xml_double_or_default("ecm_production_rate_base", 0.0);
    const double ecm_production_rate_boosted =
        read_xml_double_or_default("ecm_production_rate_boosted", 0.0);
    const double mmp2_degradation_rate =
        read_xml_double_or_default("mmp2_degradation_rate", 0.0);
    const double ecm_natural_decay_rate =
        read_xml_double_or_default("ecm_natural_decay_rate", 0.0);
    (void)ecm_natural_decay_rate;

    // PRODUCTION (stromal cells only)
    if (is_stromal && acta2_active == 1.0)
    {
        double production = ecm_production_rate_base;
        if (gli1_active == 1.0)
        {
            production = ecm_production_rate_boosted;
        }

        const double new_ecm = production * dt;
        const double new_ha  = new_ecm * ecm_ha_fraction_default;
        const double new_col = new_ecm * (1.0 - ecm_ha_fraction_default);

        const double old_density = ecm_density;
        const double old_ha_total = old_density * ecm_ha_fraction;
        const double old_col_total = old_density * (1.0 - ecm_ha_fraction);

        const double total_ha  = old_ha_total + new_ha;
        const double total_col = old_col_total + new_col;
        double new_density = total_ha + total_col;

        new_density = std::min(1.0, std::max(0.0, new_density));

        if (new_density > 0.0)
        {
            ecm_ha_fraction = total_ha / new_density;
        }
        ecm_density = new_density;
    }

    // DEGRADATION BY MMP2 (tumor cells only)
    if (is_tumor && mmp2_active == 1.0)
    {
        const double degradation = mmp2_degradation_rate * dt;
        ecm_density = ecm_density - degradation;
    }

    // TODO: Natural ECM decay should be applied once per voxel per step.
    // For now, we skip per-cell natural decay and handle turnover in ECM PDE decay.

    // EXTERNAL INTERVENTIONS
    if (intervention_state.ha_degrade_active)
    {
        const double ha_loss = intervention_state.ha_degrade_strength * dt;
        const double old_ha = ecm_density * ecm_ha_fraction;
        const double new_ha_val = std::max(0.0, old_ha - ha_loss);
        const double old_col_val = ecm_density * (1.0 - ecm_ha_fraction);
        ecm_density = new_ha_val + old_col_val;
        if (ecm_density > 0.0)
        {
            ecm_ha_fraction = new_ha_val / ecm_density;
        }
    }

    if (intervention_state.col_degrade_active)
    {
        const double col_loss = intervention_state.col_degrade_strength * dt;
        const double old_ha_val2 = ecm_density * ecm_ha_fraction;
        const double old_col_val2 = ecm_density * (1.0 - ecm_ha_fraction);
        const double new_col_val = std::max(0.0, old_col_val2 - col_loss);
        ecm_density = old_ha_val2 + new_col_val;
        if (ecm_density > 0.0)
        {
            ecm_ha_fraction = old_ha_val2 / ecm_density;
        }
    }

    ecm_density = std::max(0.0, ecm_density);
    ecm_density = std::min(1.0, ecm_density);
    ecm_ha_fraction = std::max(0.0, ecm_ha_fraction);
    ecm_ha_fraction = std::min(1.0, ecm_ha_fraction);

    rho[ecm_index] = ecm_density;
    g_ecm_ha_fraction_by_voxel[voxel_index] = ecm_ha_fraction;
}

void module8_mechanical_compaction(Cell* pCell, Phenotype& phenotype, double dt, ModulePhase phase)
{
    (void)phenotype;
    assert_expected_phase("module8_mechanical_compaction", phase, ModulePhase::WRITE);
    assert_microenvironment_write_allowed(phase, "module8_mechanical_compaction");
    log_module_trace("module8_mechanical_compaction", phase, {"ecm_density", "ecm_ha_fraction"}, {"ecm_density"});

    if (pCell == NULL) return;
    if (ecm_index < 0) return;

    const int voxel_index = get_cell_voxel_index(pCell);
    if (voxel_index < 0 ||
        voxel_index >= static_cast<int>(microenvironment.number_of_voxels()))
    {
        return;
    }

    std::vector<double>& rho = microenvironment.density_vector(voxel_index);
    if (ecm_index >= static_cast<int>(rho.size())) return;

    const double cell_density_proxy = pCell->state.simple_pressure;
    const double ecm_at_voxel = rho[ecm_index];
    const double collagen_fraction = 1.0 - get_ecm_ha_fraction(voxel_index);
    const double ecm_stiffness = ecm_at_voxel * collagen_fraction;
    const double mechanical_compaction_strength =
        read_xml_double_or_default("mechanical_compaction_strength", 0.0);
    const double solid_stress =
        cell_density_proxy * ecm_stiffness * mechanical_compaction_strength;

    set_custom_data_if_present(pCell, "mechanical_pressure", solid_stress);

    if (solid_stress > 0.0 && ecm_at_voxel > 0.0)
    {
        const double compaction_ecm_increment =
            read_xml_double_or_default("compaction_ecm_increment", 0.0);
        const double compaction = solid_stress * compaction_ecm_increment * dt;
        rho[ecm_index] = std::min(1.0, ecm_at_voxel + compaction);
    }
}

void ecm_dependent_diffusion(double dt)
{
    update_ecm_effective_diffusion_coefficients(microenvironment);
    apply_ecm_dependent_modifiers(microenvironment, dt);
}

void ecm_dependent_diffusion_solver(Microenvironment& M, double dt)
{
    // Recompute voxel-wise effective diffusion coefficients after cell modules
    // have updated ECM, then run the PDE solver and apply local barrier scaling.
    update_ecm_effective_diffusion_coefficients(M);

    if (g_base_diffusion_decay_solver != NULL &&
        g_base_diffusion_decay_solver != ecm_dependent_diffusion_solver)
    {
        g_base_diffusion_decay_solver(M, dt);
    }
    else
    {
        // Fallback: pick a default solver if none was set for some reason.
        M.auto_choose_diffusion_decay_solver();
        if (M.diffusion_decay_solver != NULL &&
            M.diffusion_decay_solver != ecm_dependent_diffusion_solver)
        {
            g_base_diffusion_decay_solver = M.diffusion_decay_solver;
            g_base_diffusion_decay_solver(M, dt);
        }
    }

    // Apply ECM-driven barrier modifiers voxel-wise.
    apply_ecm_dependent_modifiers(M, dt);
}

void register_ecm_dependent_diffusion_solver(void)
{
    if (microenvironment.diffusion_decay_solver == NULL)
    {
        microenvironment.auto_choose_diffusion_decay_solver();
    }

    if (microenvironment.diffusion_decay_solver == ecm_dependent_diffusion_solver)
    {
        return;
    }

    g_base_diffusion_decay_solver = microenvironment.diffusion_decay_solver;
    microenvironment.diffusion_decay_solver = ecm_dependent_diffusion_solver;
}

void setup_microenvironment(void)
{
    const std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();

    initialize_microenvironment();

    oxygen_index = microenvironment.find_density_index("oxygen");
    tgfb_index   = microenvironment.find_density_index("tgfb");
    shh_index    = microenvironment.find_density_index("shh");
    ecm_index    = microenvironment.find_density_index("ecm_density");
    drug_index   = microenvironment.find_density_index("drug");

    reset_ecm_ha_fraction_field(read_xml_double_or_default("ecm_ha_fraction_default", 0.5));

    register_ecm_dependent_diffusion_solver();

    const ThresholdConfig cfg = ThresholdConfig::load_from_xml();
    set_threshold_config(cfg);

    try
    {
        validate_knob_partition();
        load_knob_profiles_if_needed();

        std::string profile_name = "AsPC-1";
        if (parameters.strings.find_index("calibration_profile") >= 0)
        {
            profile_name = parameters.strings("calibration_profile");
        }

        TumorCalibrationKnobs knobs = select_profile(g_knob_profiles, profile_name);
        apply_xml_knob_overrides(knobs);
        set_active_calibration_knobs(knobs);
        std::cerr << "[setup_microenvironment] Loaded calibration profile: "
                  << profile_name << "\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "[setup_microenvironment] FATAL: failed to initialize calibration knobs: "
                  << e.what() << "\n";
        exit(EXIT_FAILURE);
    }

    const std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    std::cerr << "[setup_microenvironment] Substrate indices: "
              << "oxygen=" << oxygen_index << " tgfb=" << tgfb_index
              << " shh=" << shh_index << " ecm_density=" << ecm_index
              << " drug=" << drug_index << "\n";
    std::cerr << "[setup_microenvironment] Completed in " << ms << " ms\n";
}

void create_cell_types(void)
{
    const std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();

    if (parameters.ints.find_index("random_seed") >= 0)
    {
        SeedRandom(parameters.ints("random_seed"));
    }

    initialize_default_cell_definition();
    cell_defaults.phenotype.secretion.sync_to_microenvironment(&microenvironment);

    cell_defaults.functions.volume_update_function = standard_volume_update_function;
    cell_defaults.functions.update_velocity = standard_update_cell_velocity;
    cell_defaults.functions.update_migration_bias = NULL;

    // Register dispatcher as requested.
    cell_defaults.functions.update_phenotype = custom_function;
    cell_defaults.functions.custom_cell_rule = NULL;
    cell_defaults.functions.contact_function = NULL;
    cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
    cell_defaults.functions.calculate_distance_to_membrane = NULL;

    initialize_cell_definitions_from_pugixml();
    build_cell_definitions_maps();
    setup_signal_behavior_dictionaries();
    setup_cell_rules();

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    Cell_Definition* pStroma = find_cell_definition("stromal_cell");

    if (pTumor != NULL)
    {
        pTumor->functions.update_phenotype = custom_function;
        pTumor->functions.cell_division_function = tumor_division_callback;
        g_dispatch_by_type[pTumor->type] = tumor_phenotype_update;

        // Actor 1 identity.
        ensure_custom_scalar(pTumor, "cell_id", "none", -1.0);
        ensure_custom_scalar(pTumor, "parent_id", "none", -1.0);

        // Actor 1 genotype state (categorical: 0=WT, 1=GOF, 2=LOF).
        ensure_custom_scalar(pTumor, "gene_state_KRAS", "none", 1.0);
        ensure_custom_scalar(pTumor, "gene_state_TP53", "none", 2.0);
        ensure_custom_scalar(pTumor, "gene_state_CDKN2A", "none", 2.0);
        ensure_custom_scalar(pTumor, "gene_state_SMAD4", "none", 2.0);
        ensure_custom_scalar(pTumor, "gene_state_BCL_XL", "none", 1.0);
        ensure_custom_scalar(pTumor, "gene_state_MYC", "none", 0.0);
        ensure_custom_scalar(pTumor, "gene_state_EGFR", "none", 0.0);
        ensure_custom_scalar(pTumor, "gene_state_RB1", "none", 0.0);
        ensure_custom_scalar(pTumor, "gene_state_BAX", "none", 0.0);
        ensure_custom_scalar(pTumor, "gene_state_ZEB1", "none", 0.0);
        ensure_custom_scalar(pTumor, "gene_state_CDH1", "none", 0.0);
        ensure_custom_scalar(pTumor, "gene_state_MMP2", "none", 0.0);
        ensure_custom_scalar(pTumor, "gene_state_TGFB1", "none", 0.0);
        ensure_custom_scalar(pTumor, "gene_state_SHH", "none", 0.0);
        ensure_custom_scalar(pTumor, "gene_state_HIF1A", "none", 0.0);
        ensure_custom_scalar(pTumor, "gene_state_NRF2", "none", 0.0);
        ensure_custom_scalar(pTumor, "gene_state_ABCB1", "none", 0.0);
        ensure_custom_scalar(pTumor, "gene_state_HAS2", "none", 0.0);
        ensure_custom_scalar(pTumor, "gene_state_COL1A1", "none", 0.0);
        ensure_custom_scalar(pTumor, "gene_state_GLI1", "none", 0.0);
        ensure_custom_scalar(pTumor, "gene_state_ACTA2", "none", 0.0);

        // Ensure custom outputs used by tumor phenotype code exist.
        ensure_custom_scalar(pTumor, "is_mesenchymal", "dimensionless", 0.0);
        ensure_custom_scalar(pTumor, "drug_sensitivity", "dimensionless", 1.0);
        ensure_custom_scalar(pTumor, "hif1a_active", "dimensionless", 0.0);
        ensure_custom_scalar(pTumor, "zeb1_active", "dimensionless", 0.0);
        ensure_custom_scalar(pTumor, "cdh1_expressed", "dimensionless", 1.0);
        ensure_custom_scalar(pTumor, "mmp2_active", "dimensionless", 0.0);
        ensure_custom_scalar(pTumor, "emt_extent", "dimensionless", 0.0);
        ensure_custom_scalar(pTumor, "abcb1_active", "dimensionless", 0.0);
        ensure_custom_scalar(pTumor, "nrf2_active", "dimensionless", 0.0);
        ensure_custom_scalar(pTumor, "intracellular_drug", "dimensionless", 0.0);
        ensure_custom_scalar(pTumor, "mechanical_pressure", "dimensionless", 0.0);
        ensure_custom_scalar(pTumor, "time_since_drug_exposure", "hour", -1.0);
        ensure_custom_scalar(pTumor, "time_alive", "min", 0.0);
    }
    else
    {
        std::cerr << "[create_cell_types] WARNING: 'tumor_cell' not found in XML definitions.\n";
    }

    if (pStroma != NULL)
    {
        pStroma->functions.update_phenotype = custom_function;
        pStroma->functions.cell_division_function = stromal_division_callback;
        g_dispatch_by_type[pStroma->type] = stromal_phenotype_update;

        // Replace inherited tumor-centric schema with Actor 2 stromal schema.
        rebuild_stromal_custom_schema(pStroma, pTumor);

        // Ensure required Actor 2 stromal fields exist.
        ensure_custom_scalar(pStroma, "cell_id", "none", -1.0);
        ensure_custom_scalar(pStroma, "parent_id", "none", -1.0);
        ensure_custom_scalar(pStroma, "gene_state_SMAD4", "none", 0.0);
        ensure_custom_scalar(pStroma, "gene_state_ACTA2", "none", 0.0);
        ensure_custom_scalar(pStroma, "gene_state_GLI1", "none", 0.0);
        ensure_custom_scalar(pStroma, "gene_state_HAS2", "none", 0.0);
        ensure_custom_scalar(pStroma, "gene_state_COL1A1", "none", 0.0);
        ensure_custom_scalar(pStroma, "gene_state_TGFB1", "none", 0.0);
        ensure_custom_scalar(pStroma, "gene_state_SHH", "none", 0.0);
        ensure_custom_scalar(pStroma, "hif1a_active", "dimensionless", 0.0);
        ensure_custom_scalar(pStroma, "acta2_active", "dimensionless", 0.0);
        ensure_custom_scalar(pStroma, "gli1_active", "dimensionless", 0.0);
        ensure_custom_scalar(pStroma, "activation_mode", "dimensionless", 0.0);
        ensure_custom_scalar(pStroma, "tgfb_secretion_active", "dimensionless", 0.0);
        ensure_custom_scalar(pStroma, "ecm_production_rate", "dimensionless", 0.0);
        ensure_custom_scalar(pStroma, "time_alive", "min", 0.0);
        ensure_custom_scalar(pStroma, "mechanical_pressure", "dimensionless", 0.0);
    }
    else
    {
        std::cerr << "[create_cell_types] WARNING: 'stromal_cell' not found in XML definitions.\n";
    }

    display_cell_definitions(std::cout);

    const std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cerr << "[create_cell_types] Completed in " << ms << " ms\n";
}

void setup_tissue(void)
{
    const std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();

    // §8.1  Run baseline behavioral validation before placing any cells.
    // Failures indicate a broken gene-network invariant; do not proceed.
    if (!run_baseline_validation())
    {
        std::cerr << "[setup_tissue] FATAL: baseline validation failed."
                     "  Fix gene-network before running simulation.\n";
        exit(EXIT_FAILURE);
    }

    Cell_Definition* pTumor = find_cell_definition("tumor_cell");
    Cell_Definition* pStroma = find_cell_definition("stromal_cell");

    if (pTumor == NULL || pStroma == NULL)
    {
        std::cerr << "[setup_tissue] ERROR: tumor_cell or stromal_cell definition missing.\n";
        return;
    }

    int n_tumor = 50;
    int n_stroma = 200;
    if (parameters.ints.find_index("number_of_tumor_cells") >= 0)
    {
        n_tumor = parameters.ints("number_of_tumor_cells");
    }
    if (parameters.ints.find_index("number_of_stromal_cells") >= 0)
    {
        n_stroma = parameters.ints("number_of_stromal_cells");
    }

    // Requested default cluster scale (~50 micron). XML can override when provided.
    double tumor_cluster_radius = 50.0;
    if (parameters.doubles.find_index("tumor_cluster_radius") >= 0)
    {
        tumor_cluster_radius = parameters.doubles("tumor_cluster_radius");
    }
    if (tumor_cluster_radius <= 0.0) tumor_cluster_radius = 50.0;

    const double Xmin = microenvironment.mesh.bounding_box[0];
    const double Ymin = microenvironment.mesh.bounding_box[1];
    const double Zmin = microenvironment.mesh.bounding_box[2];
    const double Xmax = microenvironment.mesh.bounding_box[3];
    const double Ymax = microenvironment.mesh.bounding_box[4];
    const double Zmax = microenvironment.mesh.bounding_box[5];

    const double x_center = 0.5 * (Xmin + Xmax);
    const double y_center = 0.5 * (Ymin + Ymax);
    const double z_center = default_microenvironment_options.simulate_2D ? 0.0 : 0.5 * (Zmin + Zmax);

    const double x_range = Xmax - Xmin;
    const double y_range = Ymax - Ymin;
    const double z_range = default_microenvironment_options.simulate_2D ? 0.0 : (Zmax - Zmin);

    const double two_pi = 6.28318530717958647692;

    // ---- Place tumor cells in central circular cluster ----
    for (int i = 0; i < n_tumor; ++i)
    {
        const double r = tumor_cluster_radius * std::sqrt(UniformRandom());
        const double theta = two_pi * UniformRandom();

        std::vector<double> position(3, 0.0);
        position[0] = x_center + r * std::cos(theta);
        position[1] = y_center + r * std::sin(theta);
        position[2] = z_center;

        Cell* pCell = create_cell(*pTumor);
        pCell->assign_position(position);
        assign_tumor_identity_if_unset(pCell);
        initialize_boolean_network_for_cell(pCell, CellType::TUMOR);
    }

    // ---- Place stromal cells uniformly throughout the domain ----
    for (int i = 0; i < n_stroma; ++i)
    {
        std::vector<double> position(3, 0.0);
        position[0] = Xmin + UniformRandom() * x_range;
        position[1] = Ymin + UniformRandom() * y_range;
        position[2] = default_microenvironment_options.simulate_2D
                    ? 0.0
                    : (Zmin + UniformRandom() * z_range);

        Cell* pCell = create_cell(*pStroma);
        pCell->assign_position(position);
        assign_stromal_identity_if_unset(pCell);
    }

    load_interventions_at_setup();

    const std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cerr << "[setup_tissue] Placed tumor=" << n_tumor
              << " stromal=" << n_stroma
              << " cluster_radius=" << tumor_cluster_radius
              << "um in " << ms << " ms\n";
}

void custom_function(Cell* pCell, Phenotype& phenotype, double dt)
{
    if (pCell == NULL) return;
    if (dt <= 0.0) return;

    Cell_Definition* pTumorDef = find_cell_definition("tumor_cell");
    Cell_Definition* pStromaDef = find_cell_definition("stromal_cell");
    if (pTumorDef != NULL && pCell->type == pTumorDef->type)
    {
        assign_tumor_identity_if_unset(pCell);
    }
    if (pStromaDef != NULL && pCell->type == pStromaDef->type)
    {
        assign_stromal_identity_if_unset(pCell);
        // Rule 21: stromal SMAD4 is immutable WT.
        set_custom_data_if_present(pCell, "gene_state_SMAD4", 0.0);
    }

    if (phenotype.death.dead)
    {
        phenotype.secretion.set_all_secretion_to_zero();
        phenotype.secretion.set_all_uptake_to_zero();
        set_custom_data_if_present(pCell, "mechanical_pressure", 0.0);
        return;
    }

    const std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();

    // SENSING PHASE (reads environment, writes nothing)
    module1_oxygen_sensing(pCell, phenotype, dt, ModulePhase::SENSING);
    module3_stromal_activation(pCell, phenotype, dt, ModulePhase::SENSING);
    module7_drug_response(pCell, phenotype, dt, ModulePhase::SENSING);

    // DECISION PHASE (reads module outputs + environment, writes nothing to fields)
    module5_emt_engine(pCell, phenotype, dt, ModulePhase::DECISION);
    module4_proliferation_death(pCell, phenotype, dt, ModulePhase::DECISION);

    // WRITE PHASE (writes to microenvironment fields)
    module2_paracrine_secretion(pCell, phenotype, dt, ModulePhase::WRITE);
    module6_ecm_production(pCell, phenotype, dt, ModulePhase::WRITE);
    module8_mechanical_compaction(pCell, phenotype, dt, ModulePhase::WRITE);

    const std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    const long long elapsed_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

    // Thread-safe lightweight profiling counters.
    static std::atomic<unsigned long long> call_count(0ULL);
    static std::atomic<long long> cumulative_ns(0LL);

    const unsigned long long current_count =
        call_count.fetch_add(1ULL, std::memory_order_relaxed) + 1ULL;
    cumulative_ns.fetch_add(elapsed_ns, std::memory_order_relaxed);

    // Periodic lightweight profiling.
    if (current_count % 50000ULL == 0ULL)
    {
        const long long total_ns = cumulative_ns.load(std::memory_order_relaxed);
        const double avg_us =
            static_cast<double>(total_ns) / static_cast<double>(current_count) / 1000.0;

        std::cerr << "[perf] custom_function avg = "
                  << avg_us
                  << " us over " << current_count << " calls\n";
    }
}

std::vector<std::string> my_coloring_function(Cell* pCell)
{
    return paint_by_number_cell_coloring(pCell);
}
