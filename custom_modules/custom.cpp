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
#include <random>
#include <sstream>
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
std::vector<double> g_pre_diffusion_drug_by_voxel;
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

inline double relax_toward(double current, double target, double dt, double tau)
{
    const double safe_tau = std::max(1e-12, tau);
    const double fraction = clamp_unit(dt / safe_tau);
    return clamp_unit(current + (target - current) * fraction);
}

inline double smooth_threshold_response(double signal, double threshold)
{
    const double safe_signal = clamp_nonnegative(signal);
    const double safe_threshold = std::max(1e-12, clamp_nonnegative(threshold));
    if (safe_signal <= safe_threshold)
    {
        return 0.0;
    }

    const double excess = safe_signal - safe_threshold;
    return clamp_unit(excess / (excess + safe_threshold));
}

inline double stromal_gli1_support(double gli1_active)
{
    // ACTA2-positive CAFs persist after activation, but the strength of the
    // matrix-rich / TGF-beta-rich myCAF program should track the current GLI1
    // state. Keep untreated baseline support close to the earlier RC1-safe
    // envelope by letting even modest SHH sensing saturate the myCAF program;
    // intervention-specific matrix collapse is handled separately.
    return smooth_threshold_response(clamp_unit(gli1_active), 0.01);
}

inline double stromal_ecm_rate(double base_rate, double boosted_rate, double gli1_active)
{
    const double gli1_support = stromal_gli1_support(gli1_active);
    const double basal_fraction = 0.2;
    const double floor_rate = clamp_nonnegative(base_rate) * basal_fraction;
    const double ceiling_rate = std::max(floor_rate, clamp_nonnegative(boosted_rate));
    return floor_rate + (ceiling_rate - floor_rate) * gli1_support;
}

inline double stromal_tgfb_support(double gli1_active,
                                   bool shh_inhibition_active,
                                   double shh_inhibition_strength)
{
    double support = stromal_gli1_support(gli1_active);

    // SHH blockade should collapse matrix maintenance more than it collapses
    // all activated-CAF trophic signaling. Keep a partial TGF-beta program in
    // ACTA2+ CAFs during active RC3 intervention so barrier loss can dominate
    // without erasing the stromal niche completely.
    if (shh_inhibition_active && shh_inhibition_strength > 0.0)
    {
        const double support_floor = 0.6 * clamp_unit(shh_inhibition_strength);
        support = std::max(support, support_floor);
    }

    return clamp_unit(support);
}

// Stromal custom_data index for mechanical_pressure as laid out by
// rebuild_stromal_custom_schema (offset 16).  Name-based lookup via
// find_variable_index may resolve to the wrong slot after the schema
// rebuild, so stromal cells must use this hardcoded index.
static constexpr int STROMAL_MECH_PRESSURE_IDX = 16;

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

// Forward declaration — defined later in this translation unit.
double read_xml_double_or_default(const std::string& name, double default_value);

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

    // Rule 5: reset dynamic tumor state on daughter, except resistance-memory
    // states which are partially inherited to avoid artificial post-withdrawal
    // collapse during early regrowth divisions.
    const double resistance_inheritance_factor =
        clamp_unit(read_xml_double_or_default("resistance_inheritance_factor", 0.8));
    const double resistance_inheritance_noise_std =
        std::max(0.0, read_xml_double_or_default("resistance_inheritance_noise_std", 0.05));
    const double parent_nrf2 = read_custom_or_fallback(pParent, "nrf2_active", 0.0);
    const double parent_abcb1 = read_custom_or_fallback(pParent, "abcb1_active", 0.0);
    const double parent_time_since_exposure =
        read_custom_or_fallback(pParent, "time_since_drug_exposure", -1.0);
    const double parent_time_since_withdrawal =
        read_custom_or_fallback(pParent, "time_since_drug_withdrawal", -1.0);
    const double parent_in_resistance_memory =
        read_custom_or_fallback(pParent, "in_resistance_memory", 0.0);
    const double parent_nrf2_anchor =
        read_custom_or_fallback(pParent, "nrf2_withdrawal_anchor", 0.0);
    const double parent_abcb1_anchor =
        read_custom_or_fallback(pParent, "abcb1_withdrawal_anchor", 0.0);
    const auto inherited_with_noise = [&](double parent_value) -> double
    {
        double noisy =
            parent_value * resistance_inheritance_factor +
            NormalRandom(0.0, resistance_inheritance_noise_std);
        noisy = std::max(noisy, 0.15); // HARD FLOOR: Prevent lineage amnesia
        return clamp_unit(noisy);
    };

    write_custom_if_present(pChild, "hif1a_active", 0.0);
    write_custom_if_present(pChild, "zeb1_active", 0.0);
    write_custom_if_present(pChild, "cdh1_expressed", 1.0);
    write_custom_if_present(pChild, "nrf2_active", inherited_with_noise(parent_nrf2));
    {
        // Protein inheritance floor: daughter ABCB1 cannot round down to zero
        // when the parent is in post-withdrawal resistance memory.
        const bool parent_in_memory = (parent_in_resistance_memory > 0.5);
        const double base_inherited = inherited_with_noise(parent_abcb1);
        const double abcb1_floor = parent_in_memory
            ? clamp_unit(read_xml_double_or_default("abcb1_inheritance_floor", 0.15))
            : 0.0;
        write_custom_if_present(pChild, "abcb1_active", std::max(base_inherited, abcb1_floor));
    }
    write_custom_if_present(pChild, "mmp2_active", 0.0);
    write_custom_if_present(pChild, "intracellular_drug", 0.0);
    write_custom_if_present(
        pChild,
        "time_since_drug_exposure",
        (parent_time_since_exposure >= 0.0) ? parent_time_since_exposure : -1.0);
    write_custom_if_present(
        pChild,
        "time_since_drug_withdrawal",
        (parent_time_since_withdrawal >= 0.0) ? parent_time_since_withdrawal : -1.0);
    write_custom_if_present(pChild, "in_resistance_memory", parent_in_resistance_memory > 0.5 ? 1.0 : 0.0);
    write_custom_if_present(
        pChild,
        "nrf2_withdrawal_anchor",
        clamp_unit(parent_nrf2_anchor * resistance_inheritance_factor));
    write_custom_if_present(
        pChild,
        "abcb1_withdrawal_anchor",
        clamp_unit(parent_abcb1_anchor * resistance_inheritance_factor));
    write_custom_if_present(pChild, "emt_extent", 0.0);
    write_custom_if_present(pChild, "emt_activation_time", 0.0);
    write_custom_if_present(pChild, "emt_signal_time", 0.0);
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

    // XML flags control which substrates get the dt-aware correction.
    // dt_correct_drug_impedance_only (default 1): fix drug only, legacy for others.
    // dt_correct_all_substrates      (default 0): override — fix all substrates.
    static int s_mode = -1;
    if (s_mode < 0)
    {
        double flag_all = 0.0, flag_drug = 1.0;
        if (parameters.doubles.find_index("dt_correct_all_substrates") >= 0)
            flag_all = parameters.doubles("dt_correct_all_substrates");
        if (parameters.doubles.find_index("dt_correct_drug_impedance_only") >= 0)
            flag_drug = parameters.doubles("dt_correct_drug_impedance_only");
        if (flag_all > 0.5)
            s_mode = 2;  // dt-aware for all
        else if (flag_drug > 0.5)
            s_mode = 1;  // dt-aware for drug only (default)
        else
            s_mode = 0;  // legacy for all
    }

    auto base_diffusion_or_zero = [&M](int substrate_index) -> double
    {
        if (substrate_index < 0) return 0.0;
        if (substrate_index >= static_cast<int>(M.diffusion_coefficients.size())) return 0.0;
        return std::max(0.0, M.diffusion_coefficients[substrate_index]);
    };

    const bool dt_correct_drug = (s_mode >= 1);
    const bool dt_correct_others = (s_mode >= 2);

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

        // Oxygen, TGF-β, SHH: legacy form (RC1-calibrated) unless dt_correct_all
        if (base_o2 > 0.0 && n < g_effective_o2_diffusion_by_voxel.size())
        {
            const double ratio = clamp_unit(g_effective_o2_diffusion_by_voxel[n] / base_o2);
            rho[oxygen_index] = clamp_nonnegative(
                rho[oxygen_index] * (dt_correct_others ? std::pow(ratio, dt) : ratio));
        }
        if (base_tgfb > 0.0 && n < g_effective_tgfb_diffusion_by_voxel.size())
        {
            const double ratio = clamp_unit(g_effective_tgfb_diffusion_by_voxel[n] / base_tgfb);
            rho[tgfb_index] = clamp_nonnegative(
                rho[tgfb_index] * (dt_correct_others ? std::pow(ratio, dt) : ratio));
        }
        if (base_shh > 0.0 && n < g_effective_shh_diffusion_by_voxel.size())
        {
            const double ratio = clamp_unit(g_effective_shh_diffusion_by_voxel[n] / base_shh);
            rho[shh_index] = clamp_nonnegative(
                rho[shh_index] * (dt_correct_others ? std::pow(ratio, dt) : ratio));
        }

        // Drug: dt-aware correction (fixes exponential sink at small dt_diffusion)
        if (base_drug > 0.0 && n < g_effective_drug_diffusion_by_voxel.size())
        {
            const double ratio = clamp_unit(g_effective_drug_diffusion_by_voxel[n] / base_drug);
            if (dt_correct_drug &&
                n < g_pre_diffusion_drug_by_voxel.size())
            {
                const double transport_factor = clamp_unit(std::pow(ratio, 2.5));
                const double pre_drug = g_pre_diffusion_drug_by_voxel[n];
                rho[drug_index] = clamp_nonnegative(
                    pre_drug + (rho[drug_index] - pre_drug) * transport_factor);
            }
            else
            {
                rho[drug_index] = clamp_nonnegative(
                    rho[drug_index] * (dt_correct_drug ? std::pow(ratio, dt) : ratio));
            }
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
    double effective_local_shh = local_shh;
    const double shh_inhibition_start_time =
        read_xml_double_or_default("shh_inhibition_start_time", 1e18);
    const double shh_inhibition_strength =
        clamp_unit(read_xml_double_or_default("shh_inhibition_strength", 0.0));
    const bool shh_inhibition_active =
        (PhysiCell_globals.current_time >= shh_inhibition_start_time) &&
        (shh_inhibition_strength > 0.0);
    if (shh_inhibition_active)
    {
        effective_local_shh = effective_local_shh * (1.0 - shh_inhibition_strength);
    }
    const double gli1_signal = smooth_threshold_response(effective_local_shh, shh_threshold);

    if (is_stromal && previous_acta2 == 0.0)
    {
        const double combined_signal = local_tgfb + effective_local_shh;
        if (combined_signal > activation_threshold)
        {
            updated_acta2 = 1.0; // IRREVERSIBLE
            updated_gli1 = gli1_signal;
        }
    }

    if (updated_acta2 == 1.0)
    {
        updated_acta2 = 1.0; // acta2 stays ON permanently
        updated_gli1 = gli1_signal; // GLI1 is reversible and graded
    }

    if (previous_acta2 == 1.0 && updated_acta2 == 0.0)
    {
        std::cerr << "[module3_stromal_activation] FATAL: acta2_active attempted 1->0 transition.\n";
        assert(false && "acta2_active is irreversible and cannot transition from 1.0 to 0.0");
    }

    set_custom_data_if_present(pCell, "acta2_active", updated_acta2);
    set_custom_data_if_present(pCell, "gli1_active", updated_gli1);
    set_custom_data_if_present(pCell, "gene_state_ACTA2", updated_acta2);
    set_custom_data_if_present(pCell, "gene_state_GLI1", updated_gli1);

    if (updated_acta2 == 1.0)
    {
        const double tgfb_support =
            stromal_tgfb_support(updated_gli1, shh_inhibition_active, shh_inhibition_strength);
        set_custom_data_if_present(pCell, "activation_mode", 1.0);
        set_custom_data_if_present(pCell, "tgfb_secretion_active", tgfb_support);
        set_custom_data_if_present(
            pCell,
            "ecm_production_rate",
            stromal_ecm_rate(ecm_production_rate_base, ecm_production_rate_boosted, updated_gli1));
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
    double time_since_drug_withdrawal =
        read_custom_data_value_or_default(pCell, "time_since_drug_withdrawal", -1.0);
    double in_resistance_memory =
        read_custom_data_value_or_default(pCell, "in_resistance_memory", 0.0);
    double nrf2_withdrawal_anchor =
        read_custom_data_value_or_default(pCell, "nrf2_withdrawal_anchor", 0.0);
    double abcb1_withdrawal_anchor =
        read_custom_data_value_or_default(pCell, "abcb1_withdrawal_anchor", 0.0);

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
    const double resistance_withdrawal_tau =
        read_xml_double_or_default("resistance_withdrawal_tau", 10080.0);

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
        if (time_since_drug_withdrawal < 0.0)
        {
            time_since_drug_withdrawal = 0.0;
        }
    }

    // Step 2: Intracellular drug decay (metabolic clearance — always active).
    if (time_since_drug_exposure >= 0.0 && drug_natural_decay_rate > 0.0)
    {
        intracellular_drug = intracellular_drug * (1.0 - drug_natural_decay_rate * dt);
        intracellular_drug = std::max(0.0, intracellular_drug);
    }

    // Step 3: Stress sensing — Rate-Based Integration (NRF2/ABCB1 decoupled).
    // NRF2: fast transcriptional responder  t_1/2 ~ 20 min (mRNA/TF half-life).
    // ABCB1: slow efflux pump protein       t_1/2 ~ 27 h  (protein half-life).
    // Decoupled: ABCB1 decay is intrinsic to the protein pool, not chained to NRF2.
    const double configured_drug_end_time =
        read_xml_double_or_default("drug_end_time", 1e18);
    const double minutes_since_withdrawal =
        PhysiCell_globals.current_time - configured_drug_end_time;
    const bool post_withdrawal = (minutes_since_withdrawal >= 0.0);
    const bool in_withdrawal_memory =
        (post_withdrawal &&
         resistance_withdrawal_tau > 0.0 &&
         (time_since_drug_exposure >= 0.0 || local_drug > 0.0));
    in_resistance_memory = in_withdrawal_memory ? 1.0 : 0.0;

    {
        // A. NRF2 — transcription factor (fast: t_1/2 ~ 20 min = tau * ln2)
        const double nrf2_tau = 20.0;
        const double stress = intracellular_drug + (hif1a_active > 0.5 ? 0.5 : 0.0);
        const double nrf2_production = (stress > 1e-4) ? 0.05 : 0.0;
        nrf2_active += nrf2_production * dt;
        nrf2_active *= std::exp(-dt / nrf2_tau);
        nrf2_active = clamp_unit(nrf2_active);

        // Post-withdrawal transcriptional memory trace: the NRF2 program acquired
        // during treatment persists on the time-scale of resistance_withdrawal_tau.
        // Prevents complete NRF2 collapse before the ABCB1 protein pool has cleared.
        // Floor decays from 0.80 → 0 over ~7 days, staying > 0.5 for the full 72h
        // Stage-1 window (exp(-4320/10080) * 0.80 = 0.52).
        if (in_withdrawal_memory)
        {
            const double withdrawal_memory_fraction =
                std::exp(-minutes_since_withdrawal /
                         std::max(1e-12, resistance_withdrawal_tau));
            const double nrf2_memory_floor_scale =
                read_xml_double_or_default("nrf2_memory_floor_scale", 0.80);
            nrf2_active = std::max(nrf2_active,
                                   withdrawal_memory_fraction * nrf2_memory_floor_scale);
        }
    }

    // Step 4: Efflux activation (ABCB1) — protein with intrinsic half-life.
    {
        // B. ABCB1 — efflux pump protein (slow: t_1/2 ~ 27 h, tau = 1620 min)
        // Production requires significant NRF2 (> 0.5); decay is protein-intrinsic.
        // abcb1_production_rate is XML-tunable: lower values delay resistance onset,
        // giving the drug a longer kill window before ABCB1 shields go up.
        // 0.001 → SS=1.0, t_resist~10h (too fast: RC2-1 fails)
        // 0.0005 → SS=0.81, t_resist~26h (Goldilocks: 1-day kill window)
        const double abcb1_tau = 1620.0;
        const double abcb1_production_rate =
            read_xml_double_or_default("abcb1_production_rate", 0.0005);
        const double abcb1_production = (nrf2_active > 0.5) ? abcb1_production_rate : 0.0;
        abcb1_active += abcb1_production * dt;
        abcb1_active *= std::exp(-dt / abcb1_tau);
        abcb1_active = clamp_unit(abcb1_active);
    }

    // Step 5: Efflux pump action (works on intracellular drug regardless of external).
    if (intracellular_drug > 0.0 && abcb1_active > 0.0)
    {
        const double efflux_amount =
            efflux_strength * clamp_unit(abcb1_active) * intracellular_drug * dt;
        intracellular_drug = intracellular_drug - efflux_amount;
        intracellular_drug = std::max(0.0, intracellular_drug);
    }

    // Step 6: Advance exposure clock.
    if (time_since_drug_exposure >= 0.0)
    {
        time_since_drug_exposure = time_since_drug_exposure + dt;
        if (local_drug > 0.0)
        {
            nrf2_withdrawal_anchor = std::max(nrf2_withdrawal_anchor, nrf2_active);
            abcb1_withdrawal_anchor = std::max(abcb1_withdrawal_anchor, abcb1_active);
            time_since_drug_withdrawal = 0.0;
        }
        else
        {
            if (time_since_drug_withdrawal < 0.0)
            {
                time_since_drug_withdrawal = 0.0;
            }
            time_since_drug_withdrawal = time_since_drug_withdrawal + dt;
        }
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
    set_custom_data_if_present(pCell, "NRF2", nrf2_active);
    set_custom_data_if_present(pCell, "ABCB1", abcb1_active);
    set_custom_data_if_present(pCell, "time_since_drug_exposure", time_since_drug_exposure);
    set_custom_data_if_present(pCell, "time_since_drug_withdrawal", time_since_drug_withdrawal);
    set_custom_data_if_present(pCell, "in_resistance_memory", in_resistance_memory);
    set_custom_data_if_present(pCell, "nrf2_withdrawal_anchor", nrf2_withdrawal_anchor);
    set_custom_data_if_present(pCell, "abcb1_withdrawal_anchor", abcb1_withdrawal_anchor);
}

// DECISION PHASE (reads module outputs + environment, writes nothing to fields)
void module5_emt_engine(Cell* pCell, Phenotype& phenotype, double dt, ModulePhase phase)
{
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

    // Read local ECM density once (used by Fixes C, D, and ECM-EMT boost)
    const double local_ecm = (ecm_index >= 0)
        ? read_density_value(densities, ecm_index) : 0.0;

    // ── ECM-contact EMT boost (smooth ramp + cap) ───────────────────
    // Tumor cells in dense ECM (peritumoral ring) get a gradual EMT push.
    // Periphery-specific by construction: ECM only exists where CAFs
    // deposited it (ECM ~0.8–1.0 in ring), while tumor core is degraded (~0).
    //
    // Fix C addition: require an activated CAF neighbor within interaction
    // distance. This prevents self-amplifying EMT cascades in ECM zones
    // where the stroma has already been consumed or is absent.
    double ecm_boost = 0.0;
    const double ecm_emt_strength =
        read_xml_double_or_default("ecm_emt_strength", 0.0);
    if (ecm_emt_strength > 0.0)
    {
        const double ecm_emt_threshold =
            read_xml_double_or_default("ecm_emt_threshold", 0.3);
        const double ecm_emt_ramp =
            read_xml_double_or_default("ecm_emt_ramp", 0.3);
        const double ecm_emt_cap =
            read_xml_double_or_default("ecm_emt_cap", 0.25);
        const double ecm_emt_require_caf =
            read_xml_double_or_default("ecm_emt_require_caf_contact", 1.0);

        bool apply_ecm_boost = true;

        // Fix C gate: only apply ECM boost if an activated CAF is nearby
        if (ecm_emt_require_caf > 0.5)
        {
            apply_ecm_boost = false;
            Cell_Definition* pStromaDef = find_cell_definition("stromal_cell");
            if (pStromaDef != NULL)
            {
                for (int n = 0; n < (int)pCell->state.neighbors.size(); n++)
                {
                    Cell* pN = pCell->state.neighbors[n];
                    if (pN == NULL || pN->phenotype.death.dead) continue;
                    if (pN->type != pStromaDef->type) continue;
                    double a2 = read_custom_data_value_or_default(pN, "acta2_active", 0.0);
                    if (a2 >= 0.5)
                    {
                        apply_ecm_boost = true;
                        break;
                    }
                }
            }
        }

        if (apply_ecm_boost)
        {
            const double ecm_excess = std::max(0.0, local_ecm - ecm_emt_threshold);
            const double ecm_gate   = std::min(1.0, ecm_excess / ecm_emt_ramp);
            ecm_boost = std::min(ecm_emt_strength * ecm_gate, ecm_emt_cap);
            induction_signal += ecm_boost;
        }
    }

    // Reversion signal: paracrine signals + fractional ECM contribution.
    const double ecm_reversion_weight =
        read_xml_double_or_default("ecm_reversion_weight", 0.3);
    const double reversion_signal = induction_signal - (1.0 - ecm_reversion_weight) * ecm_boost;

    // ── Fix B: Hysteresis + time delay for EMT activation ───────────
    // EMT ON  requires signal > emt_induction_threshold sustained for
    //         emt_activation_delay minutes (prevents instant flipping).
    // EMT OFF requires signal < emt_off_threshold (hysteresis prevents
    //         flickering; off-threshold << on-threshold).
    const double emt_off_threshold =
        read_xml_double_or_default("emt_off_threshold", 0.15);
    const double emt_activation_delay =
        read_xml_double_or_default("emt_activation_delay", 360.0);
    const double emt_persistence_time =
        read_xml_double_or_default("emt_persistence_time", 7200.0);
    const double current_zeb1 =
        read_custom_data_value_or_default(pCell, "zeb1_active", 0.0);
    double emt_signal_time =
        read_custom_data_value_or_default(pCell, "emt_signal_time", 0.0);
    double emt_activation_time =
        read_custom_data_value_or_default(pCell, "emt_activation_time", 0.0);

    bool emt_on = false;

    if (current_zeb1 < 0.5)
    {
        // Currently epithelial: need signal above threshold for delay minutes
        if (induction_signal > emt_induction_threshold)
        {
            emt_signal_time += dt;
            if (emt_signal_time >= emt_activation_delay)
            {
                emt_on = true;
                // Record when EMT first activated
                if (emt_activation_time <= 0.0)
                {
                    emt_activation_time = PhysiCell_globals.current_time;
                }
            }
        }
        else
        {
            emt_signal_time = 0.0; // signal dropped: reset timer
        }
    }
    else
    {
        // Currently mesenchymal: revert only after persistence timer expires
        // AND reversion signal drops below threshold.
        double time_since_activation = PhysiCell_globals.current_time - emt_activation_time;
        if (time_since_activation > emt_persistence_time &&
            reversion_signal < emt_off_threshold)
        {
            emt_on = false;
            emt_signal_time = 0.0;
            emt_activation_time = 0.0; // reset for potential re-activation
        }
        else
        {
            emt_on = true;
            emt_signal_time = emt_activation_delay; // keep timer saturated
        }
    }
    set_custom_data_if_present(pCell, "emt_signal_time", emt_signal_time);
    set_custom_data_if_present(pCell, "emt_activation_time", emt_activation_time);

    if (emt_on)
    {
        set_custom_data_if_present(pCell, "zeb1_active", 1.0);
        set_custom_data_if_present(pCell, "cdh1_expressed", 0.0);
        set_custom_data_if_present(pCell, "ZEB1", 1.0);
        set_custom_data_if_present(pCell, "CDH1", 0.0);

        double emt_motility = motility_epithelial;
        double emt_adhesion = adhesion_epithelial;

        if (emt_phenotype_extent == 1)
        {
            emt_motility = motility_mesenchymal_low;
            emt_adhesion = adhesion_mesenchymal_low;
            set_custom_data_if_present(pCell, "mmp2_active", 0.0);
            set_custom_data_if_present(pCell, "MMP2", 0.0);
            set_custom_data_if_present(pCell, "emt_extent", 1.0);
        }

        if (emt_phenotype_extent == 2)
        {
            emt_motility = motility_mesenchymal_med;
            emt_adhesion = adhesion_mesenchymal_med;
            set_custom_data_if_present(pCell, "mmp2_active", 1.0);
            set_custom_data_if_present(pCell, "MMP2", 1.0);
            set_custom_data_if_present(pCell, "emt_extent", 2.0);
        }

        if (emt_phenotype_extent == 3)
        {
            emt_motility = motility_mesenchymal_high;
            emt_adhesion = adhesion_mesenchymal_high;
            set_custom_data_if_present(pCell, "mmp2_active", 1.0);
            set_custom_data_if_present(pCell, "MMP2", 1.0);
            set_custom_data_if_present(pCell, "emt_extent", 3.0);
        }

        // Fix D: Dense ECM impedes EMT cell motility.
        // Even mesenchymal cells cannot move freely through very dense
        // collagen matrix. Motility drops linearly to zero at the block
        // threshold (e.g. ECM >= 0.8 → motility = 0).
        const double ecm_motility_block =
            read_xml_double_or_default("ecm_motility_block_threshold", 0.8);
        if (ecm_motility_block > 0.0 && local_ecm > 0.0)
        {
            double ecm_motility_factor = std::max(0.0, 1.0 - local_ecm / ecm_motility_block);
            emt_motility *= ecm_motility_factor;
        }

        phenotype.motility.migration_speed = emt_motility;
        phenotype.mechanics.cell_cell_adhesion_strength = emt_adhesion;
    }
    else
    {
        set_custom_data_if_present(pCell, "zeb1_active", 0.0);
        set_custom_data_if_present(pCell, "cdh1_expressed", 1.0);
        set_custom_data_if_present(pCell, "mmp2_active", 0.0);
        set_custom_data_if_present(pCell, "ZEB1", 0.0);
        set_custom_data_if_present(pCell, "CDH1", 1.0);
        set_custom_data_if_present(pCell, "MMP2", 0.0);
        set_custom_data_if_present(pCell, "emt_extent", 0.0);
        phenotype.motility.migration_speed = motility_epithelial;
        phenotype.mechanics.cell_cell_adhesion_strength = adhesion_epithelial;
    }

    // Sync is_mesenchymal to mirror zeb1_active so evaluators can read either name.
    pCell->custom_data["is_mesenchymal"] = pCell->custom_data["zeb1_active"];
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
    double mechanical_pressure =
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

    // Fix: stromal schema has mechanical_pressure at hardcoded offset 16.
    // Name lookup may resolve to the wrong index after rebuild_stromal_custom_schema.
    if (is_stromal &&
        STROMAL_MECH_PRESSURE_IDX < static_cast<int>(pCell->custom_data.variables.size()))
    {
        mechanical_pressure = pCell->custom_data[STROMAL_MECH_PRESSURE_IDX];
    }

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
        set_custom_data_if_present(pCell, "final_prolif_rate", proliferation);

        const double efflux_level = clamp_unit(abcb1_active);
        const double effective_drug =
            intracellular_drug * (1.0 - efflux_drug_reduction * efflux_level);

        const double base_death_rate = (1.0 - apoptosis_resistance) * 0.001;

        double death_resistance_modifier = 0.0;
        if (hif1a_active == 1.0)
        {
            death_resistance_modifier = hypoxia_death_resistance_bonus;
        }

        const double drug_kill_rate = effective_drug * drug_kill_coefficient;

        // Fix A: EMT cells pay a death-rate cost (anoikis-like apoptosis)
        const double emt_death_increase =
            read_xml_double_or_default("emt_death_increase", 0.0);
        double emt_death_penalty = 0.0;
        if (zeb1_active == 1.0)
        {
            emt_death_penalty = emt_death_increase;
        }

        const double total_apoptosis_rate = std::max(
            0.0,
            base_death_rate * (1.0 - death_resistance_modifier) + drug_kill_rate + emt_death_penalty);

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
        // Stromal cells have intact checkpoints (SMAD4=WT, no KRAS/TP53/CDKN2A
        // mutations), so they respond to contact inhibition at a LOWER threshold
        // than tumor cells whose checkpoint machinery is broken.
        const double stromal_ci_threshold =
            read_xml_double_or_default("stromal_contact_inhibition_threshold",
                read_xml_double_or_default("contact_inhibition_threshold", 0.0));

        if (acta2_active == 1.0)
        {
            double caf_rate = caf_proliferation_rate;
            const double gli1_fraction = clamp_unit(gli1_active);
            caf_rate = caf_rate * (1.0 + (gli1_proliferation_boost - 1.0) * gli1_fraction);
            // Rule 37: CAFs are physical cells subject to contact inhibition.
            // Crowded CAFs also undergo pressure-induced apoptosis.
            if (mechanical_pressure > stromal_ci_threshold)
            {
                caf_rate = 0.0;
                phenotype.death.rates[apoptosis_index] =
                    read_xml_double_or_default("caf_crowded_death_rate", 0.0);
            }
            else
            {
                phenotype.death.rates[apoptosis_index] =
                    read_xml_double_or_default("caf_base_death_rate", 0.0);
            }
            phenotype.cycle.data.transition_rate(0, 0) = caf_rate;
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
        set_custom_data_if_present(pCell, "TGFB1_expr", clamp_unit(tgfb_out));
        set_custom_data_if_present(pCell, "SHH_expr", clamp_unit(shh_out));
    }

    if (is_stromal && acta2_active == 1.0)
    {
        const double tgfb_support =
            clamp_unit(read_custom_data_value_or_default(pCell, "tgfb_secretion_active", 0.0));
        if (tgfb_index >= 0 &&
            tgfb_index < static_cast<int>(phenotype.secretion.secretion_rates.size()))
        {
            phenotype.secretion.secretion_rates[tgfb_index] =
                caf_tgfb_secretion_rate * tgfb_support;
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
    const double stored_ecm_production_rate =
        read_custom_data_value_or_default(pCell, "ecm_production_rate", 0.0);
    const double mmp2_degradation_rate =
        read_xml_double_or_default("mmp2_degradation_rate", 0.0);
    const double ecm_natural_decay_rate =
        read_xml_double_or_default("ecm_natural_decay_rate", 0.0);
    (void)ecm_natural_decay_rate;

    // PRODUCTION (stromal cells only)
    if (is_stromal && acta2_active == 1.0)
    {
        double production = stored_ecm_production_rate;
        if (production <= 0.0)
        {
            production = stromal_ecm_rate(
                ecm_production_rate_base,
                ecm_production_rate_boosted,
                gli1_active);
        }

        const double shh_inhibition_start_time =
            read_xml_double_or_default("shh_inhibition_start_time", 1e18);
        const double shh_inhibition_strength =
            clamp_unit(read_xml_double_or_default("shh_inhibition_strength", 0.0));
        if (PhysiCell_globals.current_time >= shh_inhibition_start_time &&
            shh_inhibition_strength > 0.0)
        {
            const double gli1_support = stromal_gli1_support(gli1_active);
            const double matrix_maintenance =
                clamp_unit(1.0 - shh_inhibition_strength * (1.0 - gli1_support));
            production *= matrix_maintenance;
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

    if (is_stromal && acta2_active == 1.0)
    {
        const double shh_inhibition_start_time =
            read_xml_double_or_default("shh_inhibition_start_time", 1e18);
        const double shh_inhibition_strength =
            clamp_unit(read_xml_double_or_default("shh_inhibition_strength", 0.0));
        if (PhysiCell_globals.current_time >= shh_inhibition_start_time &&
            shh_inhibition_strength > 0.0)
        {
            const double gli1_support = stromal_gli1_support(gli1_active);
            const double maintenance_loss =
                ecm_natural_decay_rate * shh_inhibition_strength * (1.0 - gli1_support) * dt;
            ecm_density = ecm_density - maintenance_loss;
        }
    }

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
    const double crowding_base =
        read_xml_double_or_default("crowding_base_pressure", 0.3);
    const double solid_stress =
        cell_density_proxy * (crowding_base + ecm_stiffness * mechanical_compaction_strength);

    // Write mechanical_pressure — both cell types use the same solid_stress
    // formula so they feel identical crowding physics.  Tumor uses
    // contact_inhibition_threshold (5.0), stroma uses
    // stromal_contact_inhibition_threshold (0.8).
    {
        Cell_Definition* pStromaDef = find_cell_definition("stromal_cell");
        const bool cell_is_stromal =
            (pStromaDef != NULL) && (pCell->type == pStromaDef->type);
        if (cell_is_stromal)
        {
            if (STROMAL_MECH_PRESSURE_IDX <
                static_cast<int>(pCell->custom_data.variables.size()))
            {
                pCell->custom_data[STROMAL_MECH_PRESSURE_IDX] = solid_stress;
            }
        }
        else
        {
            set_custom_data_if_present(pCell, "mechanical_pressure", solid_stress);
        }
    }

    if (solid_stress > 0.0 && ecm_at_voxel > 0.0)
    {
        const double compaction_ecm_increment =
            read_xml_double_or_default("compaction_ecm_increment", 0.0);
        const double packing_headroom = clamp_unit(1.0 - ecm_at_voxel);
        const double compaction = solid_stress * compaction_ecm_increment * dt
            * collagen_fraction * packing_headroom;
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

    g_pre_diffusion_drug_by_voxel.clear();
    if (drug_index >= 0)
    {
        g_pre_diffusion_drug_by_voxel.resize(M.number_of_voxels(), 0.0);
        for (unsigned int n = 0; n < M.number_of_voxels(); ++n)
        {
            const std::vector<double>& rho = M.density_vector(static_cast<int>(n));
            if (drug_index < static_cast<int>(rho.size()))
            {
                g_pre_diffusion_drug_by_voxel[n] = rho[drug_index];
            }
        }
    }

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
    g_pre_diffusion_drug_by_voxel.clear();
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

    // Log ECM impedance correction mode at startup.
    {
        const double flag_all  = read_xml_double_or_default("dt_correct_all_substrates", 0.0);
        const double flag_drug = read_xml_double_or_default("dt_correct_drug_impedance_only", 1.0);
        const char* mode_str = "legacy (all substrates)";
        if (flag_all > 0.5)
            mode_str = "dt-aware (ALL substrates)";
        else if (flag_drug > 0.5)
            mode_str = "dt-aware DRUG ONLY (O2/TGFb/SHH legacy)";
        std::cerr << "[setup_microenvironment] ECM impedance correction mode: "
                  << mode_str << "\n";
    }

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
        ensure_custom_scalar(pTumor, "final_prolif_rate", "1/min", 0.0);
        ensure_custom_scalar(pTumor, "time_since_drug_exposure", "hour", -1.0);
        ensure_custom_scalar(pTumor, "time_since_drug_withdrawal", "hour", -1.0);
        ensure_custom_scalar(pTumor, "in_resistance_memory", "dimensionless", 0.0);
        ensure_custom_scalar(pTumor, "nrf2_withdrawal_anchor", "dimensionless", 0.0);
        ensure_custom_scalar(pTumor, "abcb1_withdrawal_anchor", "dimensionless", 0.0);
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

    // ---- Place stromal cells in annular ring around tumor center ----
    double stroma_inner = 100.0;
    double stroma_outer = 300.0;
    if (parameters.doubles.find_index("stroma_inner_radius") >= 0)
    {
        stroma_inner = parameters.doubles("stroma_inner_radius");
    }
    if (parameters.doubles.find_index("stroma_outer_radius") >= 0)
    {
        stroma_outer = parameters.doubles("stroma_outer_radius");
    }
    if (stroma_inner < 0.0) stroma_inner = 0.0;
    if (stroma_outer <= stroma_inner) stroma_outer = stroma_inner + 200.0;

    const double r_inner2 = stroma_inner * stroma_inner;
    const double r_outer2 = stroma_outer * stroma_outer;

    for (int i = 0; i < n_stroma; ++i)
    {
        // Uniform sampling in annulus: r = sqrt(U * (R2^2 - R1^2) + R1^2)
        const double r = std::sqrt(UniformRandom() * (r_outer2 - r_inner2) + r_inner2);
        const double theta = two_pi * UniformRandom();

        std::vector<double> position(3, 0.0);
        position[0] = x_center + r * std::cos(theta);
        position[1] = y_center + r * std::sin(theta);
        position[2] = z_center;

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
              << "um stroma_ring=[" << stroma_inner << "," << stroma_outer << "]um"
              << " in " << ms << " ms\n";
}

void custom_function(Cell* pCell, Phenotype& phenotype, double dt)
{
    if (pCell == NULL) return;
    if (dt <= 0.0) return;

    // ---- HOST DEATH CEILING: periodic check (once per 360-min interval) ----
    {
        static std::atomic<double> next_ceiling_check(360.0);
        static std::atomic<bool> ceiling_triggered(false);
        const double sim_time = PhysiCell_globals.current_time;

        if (!ceiling_triggered.load(std::memory_order_relaxed) &&
            sim_time >= next_ceiling_check.load(std::memory_order_relaxed))
        {
            // Only one thread performs the check
            double expected_next = next_ceiling_check.load(std::memory_order_acquire);
            if (sim_time >= expected_next &&
                next_ceiling_check.compare_exchange_strong(expected_next,
                    expected_next + 360.0, std::memory_order_acq_rel))
            {
                const double host_death_ceiling =
                    read_xml_double_or_default("host_death_ceiling", 50000.0);
                Cell_Definition* pTD = find_cell_definition("tumor_cell");
                if (pTD != NULL)
                {
                    int live_tumor = 0;
                    for (auto* c : *all_cells)
                    {
                        if (c != NULL && !c->phenotype.death.dead &&
                            c->type == pTD->type)
                        {
                            live_tumor++;
                        }
                    }
                    if (static_cast<double>(live_tumor) > host_death_ceiling)
                    {
                        std::cerr << "\n[HOST_DEATH_CEILING] t=" << sim_time
                                  << " live_tumor=" << live_tumor
                                  << " > ceiling=" << host_death_ceiling
                                  << "  — halting simulation.\n";
                        // Create marker file
                        std::string marker_path =
                            PhysiCell_settings.folder + "/HOST_DEATH.txt";
                        std::ofstream marker(marker_path);
                        if (marker.is_open())
                        {
                            marker << "HOST_DEATH_CEILING reached at t="
                                   << sim_time << " with " << live_tumor
                                   << " live tumor cells (ceiling="
                                   << host_death_ceiling << ")\n";
                            marker.close();
                        }
                        // Force simulation to end at next loop iteration
                        PhysiCell_settings.max_time = sim_time;
                        ceiling_triggered.store(true, std::memory_order_release);
                        return;
                    }
                }
            }
        }
        if (ceiling_triggered.load(std::memory_order_relaxed)) return;
    }

    // ---- FAST SCREEN: periodic post-withdrawal metrics + early termination ----
    // Activated by user parameter fast_screen_mode=1.
    // Logs CSV every 60 min after drug_end_time with: time, tumor_count,
    // zeb1_fraction, mean_cycle_rate, emt_to_epi_transitions (proxy).
    // Early terminates if fail criteria met after drug_end_time + 2000 min.
    {
        static std::atomic<double> next_screen_check(-1.0);
        static std::atomic<bool> screen_initialized(false);
        static std::atomic<bool> screen_terminated(false);
        static std::atomic<int> prev_zeb1_count(0);

        const double fast_mode =
            read_xml_double_or_default("fast_screen_mode", 0.0);
        if (fast_mode > 0.5 &&
            !screen_terminated.load(std::memory_order_relaxed))
        {
            const double sim_time = PhysiCell_globals.current_time;
            const double drug_end =
                read_xml_double_or_default("drug_end_time", 1e18);

            // Initialize the first check time
            if (!screen_initialized.load(std::memory_order_relaxed) &&
                sim_time >= drug_end)
            {
                bool exp = false;
                if (screen_initialized.compare_exchange_strong(
                        exp, true, std::memory_order_acq_rel))
                {
                    next_screen_check.store(drug_end + 60.0,
                                            std::memory_order_release);
                }
            }

            double nsc = next_screen_check.load(std::memory_order_relaxed);
            if (nsc > 0.0 && sim_time >= nsc)
            {
                double expected_next = nsc;
                if (next_screen_check.compare_exchange_strong(
                        expected_next, expected_next + 60.0,
                        std::memory_order_acq_rel))
                {
                    // Aggregate metrics across all tumor cells
                    Cell_Definition* pTD = find_cell_definition("tumor_cell");
                    int tumor_count = 0;
                    int zeb1_count = 0;
                    double sum_cycle_rate = 0.0;

                    if (pTD != NULL)
                    {
                        for (auto* c : *all_cells)
                        {
                            if (c == NULL || c->phenotype.death.dead) continue;
                            if (c->type != pTD->type) continue;
                            tumor_count++;
                            double z = read_custom_data_value_or_default(
                                c, "zeb1_active", 0.0);
                            if (z >= 0.5) zeb1_count++;
                            sum_cycle_rate +=
                                c->phenotype.cycle.data.transition_rate(0, 0);
                        }
                    }

                    double zeb1_frac = (tumor_count > 0)
                        ? static_cast<double>(zeb1_count) / tumor_count
                        : 0.0;
                    double mean_cycle = (tumor_count > 0)
                        ? sum_cycle_rate / tumor_count : 0.0;
                    // EMT→epithelial transitions proxy: decrease in ZEB1+ count
                    int prev_z = prev_zeb1_count.load(std::memory_order_relaxed);
                    int emt_to_epi = std::max(0, prev_z - zeb1_count);
                    prev_zeb1_count.store(zeb1_count,
                                          std::memory_order_release);

                    // Write CSV line
                    std::string csv_path =
                        PhysiCell_settings.folder + "/fast_screen.csv";
                    bool write_header = false;
                    {
                        std::ifstream test(csv_path);
                        if (!test.good()) write_header = true;
                    }
                    std::ofstream csv(csv_path, std::ios::app);
                    if (csv.is_open())
                    {
                        if (write_header)
                        {
                            csv << "time,tumor_count,zeb1_count,zeb1_fraction,"
                                << "mean_cycle_rate,emt_to_epi\n";
                        }
                        csv << sim_time << "," << tumor_count << ","
                            << zeb1_count << "," << zeb1_frac << ","
                            << mean_cycle << "," << emt_to_epi << "\n";
                        csv.close();
                    }

                    // Early termination check after drug_end + 2000
                    if (sim_time >= drug_end + 2000.0)
                    {
                        // Read baseline from first CSV line (drug_end+60)
                        std::ifstream rd(csv_path);
                        std::string hdr_line, first_line;
                        if (std::getline(rd, hdr_line) &&
                            std::getline(rd, first_line))
                        {
                            // Parse first data line for baseline zeb1_frac
                            // and cycle rate
                            double base_zeb1_frac = 0.0;
                            double base_cycle = 0.0;
                            int base_count = 0;
                            {
                                std::istringstream ss(first_line);
                                std::string tok;
                                int col = 0;
                                while (std::getline(ss, tok, ','))
                                {
                                    if (col == 1)
                                        base_count = std::stoi(tok);
                                    else if (col == 3)
                                        base_zeb1_frac = std::stod(tok);
                                    else if (col == 4)
                                        base_cycle = std::stod(tok);
                                    col++;
                                }
                            }

                            // FAIL: ZEB1 fraction flat/increased AND
                            //       cycle rate still suppressed
                            bool zeb1_stuck =
                                (zeb1_frac >= base_zeb1_frac - 0.01);
                            bool cycle_suppressed =
                                (base_cycle > 0.0)
                                    ? (mean_cycle <= base_cycle * 1.05)
                                    : (mean_cycle <= 0.0);
                            bool count_declining =
                                (tumor_count <= base_count);

                            if (zeb1_stuck && cycle_suppressed &&
                                count_declining)
                            {
                                std::cerr
                                    << "\n[FAST_SCREEN] t=" << sim_time
                                    << " FAIL: zeb1_frac=" << zeb1_frac
                                    << " (base=" << base_zeb1_frac << ")"
                                    << " cycle=" << mean_cycle
                                    << " (base=" << base_cycle << ")"
                                    << " count=" << tumor_count
                                    << " (base=" << base_count << ")"
                                    << " — early termination.\n";

                                std::string marker_path =
                                    PhysiCell_settings.folder +
                                    "/FAST_SCREEN_FAIL.txt";
                                std::ofstream marker(marker_path);
                                if (marker.is_open())
                                {
                                    marker << "FAST_SCREEN_FAIL at t="
                                           << sim_time << "\n"
                                           << "zeb1_frac=" << zeb1_frac
                                           << " base=" << base_zeb1_frac
                                           << "\nmean_cycle=" << mean_cycle
                                           << " base=" << base_cycle
                                           << "\ntumor_count=" << tumor_count
                                           << " base=" << base_count << "\n";
                                    marker.close();
                                }
                                PhysiCell_settings.max_time = sim_time;
                                screen_terminated.store(
                                    true, std::memory_order_release);
                            }
                        }
                    }
                }
            }
        }
    }

    // ---- DRUG SCHEDULING: timed Dirichlet BC on/off ----
    // drug_start_time: simulation time to enable drug boundary condition
    // drug_end_time:   simulation time to withdraw drug (set BC back to 0)
    // Only one thread per timestep executes each transition.
    {
        static std::atomic<bool> drug_on_triggered(false);
        static std::atomic<bool> drug_off_triggered(false);
        const double sim_time = PhysiCell_globals.current_time;

        const double drug_start = read_xml_double_or_default("drug_start_time", 1e18);
        const double drug_end   = read_xml_double_or_default("drug_end_time",   1e18);

        if (!drug_on_triggered.load(std::memory_order_relaxed) &&
            sim_time >= drug_start)
        {
            bool expected = false;
            if (drug_on_triggered.compare_exchange_strong(
                    expected, true, std::memory_order_acq_rel))
            {
                if (drug_index >= 0)
                {
                    const double drug_conc =
                        read_xml_double_or_default("drug_concentration", 1.0);
                    const int n_vox = (int)microenvironment.mesh.voxels.size();
                    for (int n = 0; n < n_vox; n++)
                    {
                        if (microenvironment.mesh.voxels[n].is_Dirichlet)
                        {
                            microenvironment.update_dirichlet_node(
                                n, drug_index, drug_conc);
                        }
                    }
                    std::cerr << "[DRUG_SCHEDULE] t=" << sim_time
                              << " Drug ON  concentration=" << drug_conc << "\n";
                }
            }
        }

        if (drug_on_triggered.load(std::memory_order_relaxed) &&
            !drug_off_triggered.load(std::memory_order_relaxed) &&
            sim_time >= drug_end)
        {
            bool expected = false;
            if (drug_off_triggered.compare_exchange_strong(
                    expected, true, std::memory_order_acq_rel))
            {
                if (drug_index >= 0)
                {
                    const int n_vox = (int)microenvironment.mesh.voxels.size();
                    // Zero Dirichlet BC value so boundary stops supplying drug.
                    // Interior density is NOT wiped — drug dissipates
                    // naturally via diffusion + cellular decay/uptake.
                    for (int n = 0; n < n_vox; n++)
                    {
                        if (microenvironment.mesh.voxels[n].is_Dirichlet)
                        {
                            microenvironment.update_dirichlet_node(
                                n, drug_index, 0.0);
                        }
                    }
                    std::cerr << "[DRUG_SCHEDULE] t=" << sim_time
                              << " Drug WITHDRAWN (BC→0, interior dissipates)\n";
                }
            }
        }
    }

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
        // Zero mechanical_pressure — use hardcoded index for stromal cells.
        const bool dead_is_stromal =
            (pStromaDef != NULL) && (pCell->type == pStromaDef->type);
        if (dead_is_stromal &&
            STROMAL_MECH_PRESSURE_IDX <
                static_cast<int>(pCell->custom_data.variables.size()))
        {
            pCell->custom_data[STROMAL_MECH_PRESSURE_IDX] = 0.0;
        }
        else
        {
            set_custom_data_if_present(pCell, "mechanical_pressure", 0.0);
        }
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

    // ---- ONE-TIME DIAGNOSTIC at t≈360 for PSC proliferation debugging ----
    {
        static std::atomic<bool> diag_done(false);
        const double sim_time = PhysiCell_globals.current_time;
        if (!diag_done.load(std::memory_order_relaxed) && sim_time >= 360.0)
        {
            bool expected = false;
            if (diag_done.compare_exchange_strong(expected, true,
                    std::memory_order_acq_rel))
            {
                std::cerr << "\n[DIAG_T360_STROMAL] === One-time stromal diagnostic at t="
                          << sim_time << " ===\n";

                Cell_Definition* pSD = find_cell_definition("stromal_cell");
                if (pSD != NULL)
                {
                    // Collect all live stromal cells
                    std::vector<Cell*> pscs;
                    std::vector<Cell*> cafs;
                    for (auto* c : *all_cells)
                    {
                        if (c == NULL || c->phenotype.death.dead) continue;
                        if (c->type != pSD->type) continue;
                        const double a2 = read_custom_data_value_or_default(c, "acta2_active", -999.0);
                        if (a2 < 0.5)
                            pscs.push_back(c);
                        else
                            cafs.push_back(c);
                    }
                    std::cerr << "[DIAG_T360_STROMAL] Total live stromal: "
                              << (pscs.size() + cafs.size())
                              << "  PSCs(acta2=0): " << pscs.size()
                              << "  CAFs(acta2=1): " << cafs.size() << "\n";

                    // Helper: read mechanical_pressure from stromal hardcoded index
                    auto read_stromal_mech_pressure = [](Cell* c) -> double {
                        if (STROMAL_MECH_PRESSURE_IDX <
                            static_cast<int>(c->custom_data.variables.size()))
                        {
                            return c->custom_data[STROMAL_MECH_PRESSURE_IDX];
                        }
                        return -999.0;
                    };

                    // Print 5 random PSCs
                    std::mt19937 rng(42);
                    if (pscs.size() > 5)
                    {
                        std::shuffle(pscs.begin(), pscs.end(), rng);
                    }
                    const size_t n_show = std::min(pscs.size(), (size_t)5);
                    for (size_t i = 0; i < n_show; ++i)
                    {
                        Cell* c = pscs[i];
                        const double rate = c->phenotype.cycle.data.transition_rate(0, 0);
                        const double a2 = read_custom_data_value_or_default(c, "acta2_active", -999.0);
                        const double am = read_custom_data_value_or_default(c, "activation_mode", -999.0);
                        const double mp = read_stromal_mech_pressure(c);
                        std::cerr << "[DIAG_T360_STROMAL] PSC #" << i
                                  << "  ID=" << c->ID
                                  << "  transition_rate(0,0)=" << rate
                                  << "  acta2_active=" << a2
                                  << "  activation_mode=" << am
                                  << "  mechanical_pressure[16]=" << mp << "\n";
                    }

                    // Also print 3 CAFs for comparison
                    const size_t n_caf_show = std::min(cafs.size(), (size_t)3);
                    for (size_t i = 0; i < n_caf_show; ++i)
                    {
                        Cell* c = cafs[i];
                        const double rate = c->phenotype.cycle.data.transition_rate(0, 0);
                        const double a2 = read_custom_data_value_or_default(c, "acta2_active", -999.0);
                        const double am = read_custom_data_value_or_default(c, "activation_mode", -999.0);
                        const double mp = read_stromal_mech_pressure(c);
                        std::cerr << "[DIAG_T360_STROMAL] CAF #" << i
                                  << "  ID=" << c->ID
                                  << "  transition_rate(0,0)=" << rate
                                  << "  acta2_active=" << a2
                                  << "  activation_mode=" << am
                                  << "  mechanical_pressure[16]=" << mp << "\n";
                    }

                    // Verify: 3 random stromal cells mechanical_pressure
                    // If nonzero, the hardcoded-index write in Module 8 is working.
                    std::vector<Cell*> all_stromal;
                    all_stromal.insert(all_stromal.end(), pscs.begin(), pscs.end());
                    all_stromal.insert(all_stromal.end(), cafs.begin(), cafs.end());
                    std::shuffle(all_stromal.begin(), all_stromal.end(), rng);
                    const size_t n_verify = std::min(all_stromal.size(), (size_t)3);
                    std::cerr << "[DIAG_T360_MECH_VERIFY] Checking mechanical_pressure for "
                              << n_verify << " random stromal cells:\n";
                    for (size_t i = 0; i < n_verify; ++i)
                    {
                        Cell* c = all_stromal[i];
                        const double mp = read_stromal_mech_pressure(c);
                        const double mp_name = read_custom_data_value_or_default(
                            c, "mechanical_pressure", -999.0);
                        std::cerr << "[DIAG_T360_MECH_VERIFY]   cell ID=" << c->ID
                                  << "  idx[16]=" << mp
                                  << "  name_lookup=" << mp_name
                                  << "  simple_pressure=" << c->state.simple_pressure
                                  << (mp > 0.0 ? "  OK" : "  STILL_ZERO") << "\n";
                    }
                }
                std::cerr << "[DIAG_T360_STROMAL] === end diagnostic ===\n\n";
            }
        }
    }
    // ---- END DIAGNOSTIC ----

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
