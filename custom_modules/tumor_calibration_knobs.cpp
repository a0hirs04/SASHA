#include "tumor_calibration_knobs.h"

#include "../BioFVM/json.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <fstream>
#include <stdexcept>
#include <unordered_set>

using json = nlohmann::json;

namespace
{

struct KnobSpec
{
    const char*          id;
    KnobCategory         category;
    KnobRangeType        range;
    KnobDangerDirection  danger;
};

static constexpr std::array<KnobSpec, 13> KNOB_SPECS = {{
    { "tgfb_secretion_rate",            KnobCategory::Targetable, KnobRangeType::Continuous,  KnobDangerDirection::HigherIsWorse },
    { "shh_secretion_rate",             KnobCategory::Targetable, KnobRangeType::Continuous,  KnobDangerDirection::ContextDependent },
    { "tgfb_brake_sensitivity",         KnobCategory::Fixed,      KnobRangeType::Continuous,  KnobDangerDirection::LowerIsWorse },
    { "proliferation_rate",             KnobCategory::Fixed,      KnobRangeType::Continuous,  KnobDangerDirection::HigherIsWorse },
    { "checkpoint_integrity",           KnobCategory::Fixed,      KnobRangeType::Continuous,  KnobDangerDirection::LowerIsWorse },
    { "apoptosis_resistance",           KnobCategory::Fixed,      KnobRangeType::Continuous,  KnobDangerDirection::HigherIsWorse },
    { "emt_induction_threshold",        KnobCategory::Observable, KnobRangeType::Continuous,  KnobDangerDirection::LowerIsWorse },
    { "emt_phenotype_extent",           KnobCategory::Observable, KnobRangeType::Categorical, KnobDangerDirection::HigherIsWorse },
    { "hypoxia_response_threshold",     KnobCategory::Observable, KnobRangeType::Continuous,  KnobDangerDirection::LowerIsWorse },
    { "hypoxia_phenotype_shift",        KnobCategory::Observable, KnobRangeType::Categorical, KnobDangerDirection::HigherIsWorse },
    { "efflux_induction_delay",         KnobCategory::Targetable, KnobRangeType::Continuous,  KnobDangerDirection::HigherIsWorse },
    { "efflux_strength",                KnobCategory::Targetable, KnobRangeType::Continuous,  KnobDangerDirection::HigherIsWorse },
    { "mechanical_compaction_strength", KnobCategory::Fixed,      KnobRangeType::Continuous,  KnobDangerDirection::HigherIsWorse },
}};

TumorCalibrationKnobs g_active_knobs = make_aspc1_defaults();

inline double clamp01(double x)
{
    if (x < 0.0) return 0.0;
    if (x > 1.0) return 1.0;
    return x;
}

std::string upper_copy(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return static_cast<char>(std::toupper(c));
    });
    return s;
}

const KnobSpec* find_spec(const std::string& id)
{
    for (const KnobSpec& spec : KNOB_SPECS)
    {
        if (id == spec.id) return &spec;
    }
    return nullptr;
}

void set_knob_numeric(TumorCalibrationKnobs& knobs, const std::string& id, double value)
{
    const double v = clamp01(value);
    if (id == "tgfb_secretion_rate")            { knobs.tgfb_secretion_rate = v; return; }
    if (id == "shh_secretion_rate")             { knobs.shh_secretion_rate = v; return; }
    if (id == "tgfb_brake_sensitivity")         { knobs.tgfb_brake_sensitivity = v; return; }
    if (id == "proliferation_rate")             { knobs.proliferation_rate = v; return; }
    if (id == "checkpoint_integrity")           { knobs.checkpoint_integrity = v; return; }
    if (id == "apoptosis_resistance")           { knobs.apoptosis_resistance = v; return; }
    if (id == "emt_induction_threshold")        { knobs.emt_induction_threshold = v; return; }
    if (id == "hypoxia_response_threshold")     { knobs.hypoxia_response_threshold = v; return; }
    if (id == "efflux_induction_delay")         { knobs.efflux_induction_delay = v; return; }
    if (id == "efflux_strength")                { knobs.efflux_strength = v; return; }
    if (id == "mechanical_compaction_strength") { knobs.mechanical_compaction_strength = v; return; }

    throw std::runtime_error("Unknown continuous knob: " + id);
}

void set_knob_categorical(TumorCalibrationKnobs& knobs, const std::string& id, KnobLevel level)
{
    if (id == "emt_phenotype_extent")
    {
        knobs.emt_phenotype_extent = level;
        return;
    }
    if (id == "hypoxia_phenotype_shift")
    {
        knobs.hypoxia_phenotype_shift = level;
        return;
    }

    throw std::runtime_error("Unknown categorical knob: " + id);
}

void apply_profile_json_overrides(TumorCalibrationKnobs& knobs, const json& profile_json)
{
    if (!profile_json.is_object())
    {
        throw std::runtime_error("Profile payload must be a JSON object.");
    }

    for (auto it = profile_json.begin(); it != profile_json.end(); ++it)
    {
        const std::string key = it.key();
        const KnobSpec* spec = find_spec(key);
        if (spec == nullptr)
        {
            continue; // allow extra metadata keys.
        }

        if (spec->range == KnobRangeType::Continuous)
        {
            if (!it.value().is_number())
            {
                throw std::runtime_error("Knob '" + key + "' must be numeric.");
            }
            set_knob_numeric(knobs, key, it.value().get<double>());
        }
        else
        {
            if (!it.value().is_string())
            {
                throw std::runtime_error("Knob '" + key + "' must be categorical string.");
            }
            set_knob_categorical(knobs, key, knob_level_from_string(it.value().get<std::string>()));
        }
    }
}

int count_recognized_profile_keys(const json& profile_json)
{
    if (!profile_json.is_object()) return 0;
    int count = 0;
    for (auto it = profile_json.begin(); it != profile_json.end(); ++it)
    {
        if (find_spec(it.key()) != nullptr) ++count;
    }
    return count;
}

} // namespace

TumorCalibrationKnobs make_aspc1_defaults()
{
    TumorCalibrationKnobs k;
    k.tgfb_secretion_rate = 0.8;
    k.shh_secretion_rate = 0.7;
    k.tgfb_brake_sensitivity = 0.05;
    k.proliferation_rate = 0.75;
    k.checkpoint_integrity = 0.05;
    k.apoptosis_resistance = 0.85;
    k.emt_induction_threshold = 0.3;
    k.emt_phenotype_extent = KnobLevel::HIGH;
    k.hypoxia_response_threshold = 0.35;
    k.hypoxia_phenotype_shift = KnobLevel::HIGH;
    k.efflux_induction_delay = 0.5;
    k.efflux_strength = 0.7;
    k.mechanical_compaction_strength = 0.6;
    return k;
}

TumorCalibrationKnobs make_panc1_defaults()
{
    TumorCalibrationKnobs k = make_aspc1_defaults();
    // SMAD4-intact control for Knob 2 behavior checks.
    k.tgfb_brake_sensitivity = 0.9;
    return k;
}

TumorCalibrationProfiles load_tumor_calibration_profiles(const std::string& json_path)
{
    TumorCalibrationProfiles profiles;
    profiles.aspc1 = make_aspc1_defaults();
    profiles.panc1 = make_panc1_defaults();

    std::ifstream ifs(json_path);
    if (!ifs.is_open())
    {
        throw std::runtime_error("Cannot open tumor calibration profiles JSON: " + json_path);
    }

    json root;
    ifs >> root;
    if (!root.is_object())
    {
        throw std::runtime_error("Calibration profile JSON root must be an object.");
    }

    json profiles_node;
    if (root.contains("profiles") && root["profiles"].is_object())
    {
        profiles_node = root["profiles"];
    }
    else
    {
        profiles_node = root;
    }

    if (!profiles_node.contains("AsPC-1"))
    {
        throw std::runtime_error("Calibration profile JSON must contain an AsPC-1 profile.");
    }
    if (profiles_node.contains("AsPC-1"))
    {
        if (count_recognized_profile_keys(profiles_node["AsPC-1"]) != 13)
        {
            throw std::runtime_error("AsPC-1 profile must define exactly 13 knobs.");
        }
        apply_profile_json_overrides(profiles.aspc1, profiles_node["AsPC-1"]);
    }
    if (profiles_node.contains("PANC-1"))
    {
        apply_profile_json_overrides(profiles.panc1, profiles_node["PANC-1"]);
    }

    validate_knob_partition();
    return profiles;
}

TumorCalibrationKnobs select_profile(const TumorCalibrationProfiles& profiles,
                                     const std::string& profile_name)
{
    if (profile_name.empty() || profile_name == "AsPC-1")
    {
        return profiles.aspc1;
    }
    if (profile_name == "PANC-1")
    {
        return profiles.panc1;
    }
    throw std::runtime_error("Unknown calibration profile: '" + profile_name
                             + "' (expected AsPC-1 or PANC-1)");
}

KnobLevel knob_level_from_string(const std::string& value)
{
    const std::string s = upper_copy(value);
    if (s == "LOW" || s == "L") return KnobLevel::LOW;
    if (s == "MEDIUM" || s == "M") return KnobLevel::MEDIUM;
    if (s == "HIGH" || s == "H") return KnobLevel::HIGH;
    throw std::runtime_error("Invalid categorical knob level: '" + value
                             + "' (expected LOW/MEDIUM/HIGH)");
}

const char* knob_level_to_string(KnobLevel level)
{
    switch (level)
    {
        case KnobLevel::LOW: return "LOW";
        case KnobLevel::MEDIUM: return "MEDIUM";
        case KnobLevel::HIGH: return "HIGH";
    }
    return "MEDIUM";
}

double knob_level_multiplier(KnobLevel level)
{
    switch (level)
    {
        case KnobLevel::LOW: return 0.75;
        case KnobLevel::MEDIUM: return 1.0;
        case KnobLevel::HIGH: return 1.25;
    }
    return 1.0;
}

bool is_known_knob(const std::string& knob)
{
    return find_spec(knob) != nullptr;
}

bool is_targetable_knob(const std::string& knob)
{
    const KnobSpec* spec = find_spec(knob);
    return spec != nullptr && spec->category == KnobCategory::Targetable;
}

bool is_observable_knob(const std::string& knob)
{
    const KnobSpec* spec = find_spec(knob);
    return spec != nullptr && spec->category == KnobCategory::Observable;
}

bool is_fixed_knob(const std::string& knob)
{
    const KnobSpec* spec = find_spec(knob);
    return spec != nullptr && spec->category == KnobCategory::Fixed;
}

KnobDangerDirection knob_danger_direction(const std::string& knob)
{
    const KnobSpec* spec = find_spec(knob);
    if (spec == nullptr)
    {
        throw std::runtime_error("Unknown knob in knob_danger_direction(): '" + knob + "'.");
    }
    return spec->danger;
}

void validate_knob_partition()
{
    int fixed_count = 0;
    int observable_count = 0;
    int targetable_count = 0;

    std::unordered_set<std::string> targetable_ids;

    for (const KnobSpec& spec : KNOB_SPECS)
    {
        if (spec.category == KnobCategory::Fixed) ++fixed_count;
        if (spec.category == KnobCategory::Observable) ++observable_count;
        if (spec.category == KnobCategory::Targetable)
        {
            ++targetable_count;
            targetable_ids.insert(spec.id);
        }
    }

    if (static_cast<int>(KNOB_SPECS.size()) != 13)
    {
        throw std::runtime_error("Knob partition validation failed: expected exactly 13 knobs.");
    }
    if (fixed_count != 5 || observable_count != 4 || targetable_count != 4)
    {
        throw std::runtime_error("Knob partition validation failed: expected counts Fixed=5, Observable=4, Targetable=4.");
    }

    const std::array<std::string, 4> expected_targetable = {{
        "tgfb_secretion_rate",
        "shh_secretion_rate",
        "efflux_induction_delay",
        "efflux_strength"
    }};
    for (const std::string& id : expected_targetable)
    {
        if (targetable_ids.find(id) == targetable_ids.end())
        {
            throw std::runtime_error("Knob partition validation failed: missing targetable knob '" + id + "'.");
        }
    }
}

void apply_targetable_knob_interventions(TumorCalibrationKnobs& knobs,
                                         const std::vector<KnobIntervention>& interventions)
{
    if (interventions.size() > 4)
    {
        throw std::runtime_error("Partition violation: at most 4 targetable knob interventions are allowed.");
    }

    std::unordered_set<std::string> seen_targetable;
    for (const KnobIntervention& iv : interventions)
    {
        if (!is_targetable_knob(iv.knob))
        {
            if (is_known_knob(iv.knob))
            {
                throw std::runtime_error("Partition violation: knob '" + iv.knob
                                         + "' is not targetable (Fixed/Observable are read-only)."
                                         );
            }
            throw std::runtime_error("Unknown knob in intervention: '" + iv.knob + "'.");
        }
        if (seen_targetable.find(iv.knob) != seen_targetable.end())
        {
            throw std::runtime_error("Partition violation: duplicate intervention for knob '" + iv.knob + "'.");
        }
        seen_targetable.insert(iv.knob);

        const std::string effect = upper_copy(iv.effect);
        const double s = clamp01(iv.strength);

        double current = 0.0;
        if (iv.knob == "tgfb_secretion_rate") current = knobs.tgfb_secretion_rate;
        else if (iv.knob == "shh_secretion_rate") current = knobs.shh_secretion_rate;
        else if (iv.knob == "efflux_induction_delay") current = knobs.efflux_induction_delay;
        else if (iv.knob == "efflux_strength") current = knobs.efflux_strength;
        else throw std::runtime_error("Internal error: non-continuous targetable knob '" + iv.knob + "'.");

        double updated = current;
        if (effect == "INHIBIT")
        {
            updated = current * (1.0 - s);
        }
        else if (effect == "ACTIVATE")
        {
            updated = current + s * (1.0 - current);
        }
        else
        {
            throw std::runtime_error("Unsupported knob intervention effect: '" + iv.effect
                                     + "' (expected INHIBIT or ACTIVATE)");
        }

        set_knob_numeric(knobs, iv.knob, updated);
    }
}

const TumorCalibrationKnobs& get_active_calibration_knobs()
{
    return g_active_knobs;
}

void set_active_calibration_knobs(const TumorCalibrationKnobs& knobs)
{
    g_active_knobs = knobs;
}
