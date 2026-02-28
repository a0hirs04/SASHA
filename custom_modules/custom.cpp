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

namespace
{

typedef void (*PhenotypeDispatchFn)(Cell*, Phenotype&, double);

std::unordered_map<int, PhenotypeDispatchFn> g_dispatch_by_type;

// Store whichever diffusion-decay solver was configured before we wrapped it.
void (*g_base_diffusion_decay_solver)(Microenvironment&, double) = NULL;

// Requested barrier constants.
static const double ECM_BLOCKING_FACTOR   = 0.8;
static const double DRUG_ECM_DECAY_COEFF  = 0.001; // local post-solve drug loss scaling
static const double OXYGEN_ECM_DECAY_COEFF= 0.001; // local post-solve oxygen perfusion penalty

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
    if (oxygen_index < 0 || drug_index < 0 || ecm_index < 0)
    {
        if (!warned_missing_substrates)
        {
            std::cerr << "[ecm_dependent_diffusion] WARNING: missing substrate indices "
                      << "(oxygen/drug/ecm_density). Skipping barrier modifier.\n";
            warned_missing_substrates = true;
        }
        return;
    }

    const unsigned int n_voxels = M.number_of_voxels();
    for (unsigned int n = 0; n < n_voxels; ++n)
    {
        std::vector<double>& rho = M.density_vector(static_cast<int>(n));
        if (ecm_index >= static_cast<int>(rho.size()) ||
            drug_index >= static_cast<int>(rho.size()) ||
            oxygen_index >= static_cast<int>(rho.size()))
        {
            continue;
        }

        const double ecm = clamp_unit(rho[ecm_index]);

        // Drug barrier model:
        //   D_drug(ecm) = D_base * (1 - 0.8 * ecm)
        // With uniform substrate diffusion in this BioFVM mode, approximate by
        // local post-diffusion attenuation.
        const double blocking = 1.0 - ECM_BLOCKING_FACTOR * ecm;
        (void)blocking;

        double drug_multiplier = 1.0 - ECM_BLOCKING_FACTOR * ecm * DRUG_ECM_DECAY_COEFF * dt;
        drug_multiplier = std::max(0.0, drug_multiplier);
        rho[drug_index] = clamp_nonnegative(rho[drug_index] * drug_multiplier);

        // Oxygen perfusion penalty in dense stroma:
        //   o2_decay_boost = 1 + 2 * ecm
        const double o2_decay_boost = 1.0 + 2.0 * ecm;
        double oxygen_multiplier = 1.0 - (o2_decay_boost - 1.0) * OXYGEN_ECM_DECAY_COEFF * dt;
        oxygen_multiplier = std::max(0.0, oxygen_multiplier);
        rho[oxygen_index] = clamp_nonnegative(rho[oxygen_index] * oxygen_multiplier);
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

} // namespace

// SENSING PHASE (reads environment, writes nothing)
void module1_oxygen_sensing(Cell* pCell, Phenotype& phenotype, double dt, ModulePhase phase)
{
    (void)phenotype;
    (void)dt;
    assert_expected_phase("module1_oxygen_sensing", phase, ModulePhase::SENSING);

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
    (void)phenotype;
    (void)dt;
    assert_expected_phase("module3_stromal_activation", phase, ModulePhase::SENSING);

    if (pCell == NULL) return;

    const std::vector<double>& densities = pCell->nearest_density_vector();
    const double local_tgfb = read_density_value(densities, tgfb_index);
    const double local_shh = read_density_value(densities, shh_index);

    const double activation_threshold =
        read_xml_double_or_default("tgfb_activation_threshold", 0.0);
    const double shh_threshold =
        read_xml_double_or_default("shh_activation_threshold", 0.0);

    const double previous_acta2 =
        read_custom_data_value_or_default(pCell, "acta2_active", 0.0);
    double updated_acta2 = previous_acta2;

    const bool is_stromal = (pCell->type_name == "stromal_cell");
    const bool is_caf = (pCell->type_name == "CAF");

    if (is_stromal && previous_acta2 == 0.0)
    {
        const double combined_signal = local_tgfb + local_shh;
        if (combined_signal > activation_threshold)
        {
            updated_acta2 = 1.0; // IRREVERSIBLE
            pCell->type_name = "CAF";

            if (local_shh > 0.0)
            {
                set_custom_data_if_present(pCell, "gli1_active", 1.0);
            }
        }
    }

    if (is_caf && previous_acta2 == 1.0)
    {
        updated_acta2 = 1.0; // acta2 stays ON permanently

        if (local_shh > shh_threshold)
        {
            set_custom_data_if_present(pCell, "gli1_active", 1.0);
        }
        else
        {
            set_custom_data_if_present(pCell, "gli1_active", 0.0); // GLI1 is reversible
        }
    }

    if (previous_acta2 == 1.0 && updated_acta2 == 0.0)
    {
        std::cerr << "[module3_stromal_activation] FATAL: acta2_active attempted 1->0 transition.\n";
        assert(false && "acta2_active is irreversible and cannot transition from 1.0 to 0.0");
    }

    set_custom_data_if_present(pCell, "acta2_active", updated_acta2);
}

void module7_drug_response(Cell* pCell, Phenotype& phenotype, double dt, ModulePhase phase)
{
    (void)pCell;
    (void)phenotype;
    (void)dt;
    assert_expected_phase("module7_drug_response", phase, ModulePhase::SENSING);
}

// DECISION PHASE (reads module outputs + environment, writes nothing to fields)
void module5_emt_engine(Cell* pCell, Phenotype& phenotype, double dt, ModulePhase phase)
{
    (void)pCell;
    (void)phenotype;
    (void)dt;
    assert_expected_phase("module5_emt_engine", phase, ModulePhase::DECISION);
}

void module4_proliferation_death(Cell* pCell, Phenotype& phenotype, double dt, ModulePhase phase)
{
    (void)pCell;
    (void)phenotype;
    (void)dt;
    assert_expected_phase("module4_proliferation_death", phase, ModulePhase::DECISION);
}

// WRITE PHASE (writes to microenvironment fields)
void module2_paracrine_secretion(Cell* pCell, Phenotype& phenotype, double dt, ModulePhase phase)
{
    (void)pCell;
    (void)phenotype;
    (void)dt;
    assert_expected_phase("module2_paracrine_secretion", phase, ModulePhase::WRITE);
    // Placeholder guard for future field writes in this module.
    assert_microenvironment_write_allowed(phase, "module2_paracrine_secretion");
}

void module6_ecm_production(Cell* pCell, Phenotype& phenotype, double dt, ModulePhase phase)
{
    (void)pCell;
    (void)phenotype;
    (void)dt;
    assert_expected_phase("module6_ecm_production", phase, ModulePhase::WRITE);
    // Placeholder guard for future field writes in this module.
    assert_microenvironment_write_allowed(phase, "module6_ecm_production");
}

void module8_mechanical_compaction(Cell* pCell, Phenotype& phenotype, double dt, ModulePhase phase)
{
    (void)pCell;
    (void)phenotype;
    (void)dt;
    assert_expected_phase("module8_mechanical_compaction", phase, ModulePhase::WRITE);
    // Placeholder guard for future field writes in this module.
    assert_microenvironment_write_allowed(phase, "module8_mechanical_compaction");
}

void ecm_dependent_diffusion(double dt)
{
    apply_ecm_dependent_modifiers(microenvironment, dt);
}

void ecm_dependent_diffusion_solver(Microenvironment& M, double dt)
{
    // Run the base diffusion/decay solver first.
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

    // Then apply ECM-driven barrier modifiers voxel-wise.
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
        g_dispatch_by_type[pTumor->type] = tumor_phenotype_update;

        // Ensure custom outputs used by tumor phenotype code exist.
        ensure_custom_scalar(pTumor, "is_mesenchymal", "dimensionless", 0.0);
        ensure_custom_scalar(pTumor, "drug_sensitivity", "dimensionless", 1.0);
        ensure_custom_scalar(pTumor, "hif1a_active", "dimensionless", 0.0);
    }
    else
    {
        std::cerr << "[create_cell_types] WARNING: 'tumor_cell' not found in XML definitions.\n";
    }

    if (pStroma != NULL)
    {
        pStroma->functions.update_phenotype = custom_function;
        g_dispatch_by_type[pStroma->type] = stromal_phenotype_update;

        // Ensure custom outputs used by stromal phenotype code exist.
        ensure_custom_scalar(pStroma, "is_activated", "dimensionless", 0.0);
        ensure_custom_scalar(pStroma, "ecm_production_rate", "dimensionless", 0.0);
        ensure_custom_scalar(pStroma, "local_ecm_density", "dimensionless", 0.0);
        ensure_custom_scalar(pStroma, "hif1a_active", "dimensionless", 0.0);
        ensure_custom_scalar(pStroma, "acta2_active", "dimensionless", 0.0);
        ensure_custom_scalar(pStroma, "gli1_active", "dimensionless", 0.0);
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
        initialize_boolean_network_for_cell(pCell, CellType::STROMA);
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
