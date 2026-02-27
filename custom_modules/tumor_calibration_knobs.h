#ifndef TUMOR_CALIBRATION_KNOBS_H
#define TUMOR_CALIBRATION_KNOBS_H

#include <string>
#include <vector>

// 13-knob calibration model from Tumor_Calibration_Knobs.pdf.

enum class KnobCategory
{
    Fixed,
    Observable,
    Targetable
};

enum class KnobRangeType
{
    Continuous,
    Categorical
};

enum class KnobDangerDirection
{
    HigherIsWorse,
    LowerIsWorse,
    ContextDependent
};

enum class KnobLevel
{
    LOW,
    MEDIUM,
    HIGH
};

struct KnobIntervention
{
    std::string knob;
    std::string effect;   // INHIBIT | ACTIVATE
    double      strength; // [0,1]
    std::string name;
};

struct TumorCalibrationKnobs
{
    // Trait 1
    double tgfb_secretion_rate = 0.8;   // 1a (Targetable)
    double shh_secretion_rate  = 0.7;   // 1b (Targetable)

    // Trait 2
    double tgfb_brake_sensitivity = 0.05; // 2 (Fixed)

    // Trait 3
    double proliferation_rate = 0.75;   // 3a (Fixed)
    double checkpoint_integrity = 0.05; // 3b (Fixed)

    // Trait 4
    double apoptosis_resistance = 0.85; // 4 (Fixed)

    // Trait 5
    double    emt_induction_threshold = 0.3;  // 5a (Observable)
    KnobLevel emt_phenotype_extent = KnobLevel::HIGH; // 5b (Observable)

    // Trait 6
    double    hypoxia_response_threshold = 0.35; // 6a (Observable)
    KnobLevel hypoxia_phenotype_shift = KnobLevel::HIGH; // 6b (Observable)

    // Trait 7
    double efflux_induction_delay = 0.5; // 7a (Targetable)
    double efflux_strength = 0.7;        // 7b (Targetable)

    // Trait 8
    double mechanical_compaction_strength = 0.6; // 8 (Fixed)
};

struct TumorCalibrationProfiles
{
    TumorCalibrationKnobs aspc1;
    TumorCalibrationKnobs panc1;
};

TumorCalibrationKnobs make_aspc1_defaults();
TumorCalibrationKnobs make_panc1_defaults();

// Throws std::runtime_error on parse/validation errors.
TumorCalibrationProfiles load_tumor_calibration_profiles(const std::string& json_path);
TumorCalibrationKnobs select_profile(const TumorCalibrationProfiles& profiles,
                                     const std::string& profile_name);

KnobLevel knob_level_from_string(const std::string& value);
const char* knob_level_to_string(KnobLevel level);

double knob_level_multiplier(KnobLevel level);

bool is_known_knob(const std::string& knob);
bool is_targetable_knob(const std::string& knob);
bool is_observable_knob(const std::string& knob);
bool is_fixed_knob(const std::string& knob);
KnobDangerDirection knob_danger_direction(const std::string& knob);

void validate_knob_partition();

// Mutates only targetable knobs. Throws std::runtime_error for any non-targetable override.
void apply_targetable_knob_interventions(TumorCalibrationKnobs& knobs,
                                         const std::vector<KnobIntervention>& interventions);

// Global runtime state (set once at setup; read by phenotype and network code).
const TumorCalibrationKnobs& get_active_calibration_knobs();
void set_active_calibration_knobs(const TumorCalibrationKnobs& knobs);

#endif // TUMOR_CALIBRATION_KNOBS_H
