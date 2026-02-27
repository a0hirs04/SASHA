#ifndef BOOLEAN_NETWORK_H
#define BOOLEAN_NETWORK_H

// ============================================================================
// boolean_network.h — Stroma World / PROJECT-NORTHSTAR
//
// Continuous Boolean network for the 21-gene PDAC tumor-stroma model.
//
// DESIGN PHILOSOPHY:
//   Classical Boolean networks use discrete {0,1} gene states. This model
//   uses continuous [0.0, 1.0] states so that:
//     (a) the evolutionary algorithm can evolve graded drug interventions
//         (e.g., "reduce BCL_XL to 0.3" rather than simply "block BCL_XL"),
//     (b) partial signaling states are representable (heterogeneous tumors),
//     (c) smooth bifurcations replace hard switches (more numerically stable).
//
//   Update rule (per gene i, per timestep dt):
//       gene[i] += dt * (target[i] - gene[i]) / tau[i]
//
//   where target[i] is computed from Boolean-logic rules evaluated on the
//   current gene state vector. This is equivalent to a first-order linear
//   filter with time constant tau, converging toward the Boolean attractor.
//
// TIME CONSTANTS (tau):
//   Signaling genes (rapid post-translational)  : tau =   6.0 min
//   Transcription factors (slower nuclear)      : tau =  60.0 min
//   Structural/secreted proteins (slow turnover): tau = 360.0 min
//
// MUTANT LOCKING:
//   PDAC-defining mutations lock certain genes permanently:
//     KRAS=1.0, TP53=0.0, CDKN2A=0.0, SMAD4=0.0 (tumor cells)
//   Locked genes are skipped in update(); their states cannot be changed
//   even by drug interventions unless is_mutant[] is explicitly cleared.
//
// USAGE:
//   ThresholdConfig cfg = ThresholdConfig::load_from_xml(); // once at setup
//   BooleanNetwork bn;
//   bn.initialize(CellType::TUMOR, params);
//   bn.sync_from_cell(pCell);                // read current PhysiCell custom_data
//   bn.update(dt, o2, tgfb, shh, drug, cfg); // run one timestep using XML thresholds
//   bn.sync_to_cell(pCell);                  // write updated states back
// ============================================================================

#include "gene_definitions.h"
#include "../core/PhysiCell_cell.h"

#include <string>
#include <vector>
#include <stdexcept>

// ----------------------------------------------------------------------------
// EffectType — how an intervention modifies a gene's state
//
//   INHIBIT / ACTIVATE  — reversible, drug-like: applied after every update()
//                         call. Removing the intervention restores dynamics.
//   KNOCKDOWN / OVEREXPRESS — permanent, gene-therapy-like: lock is_mutant[i]
//                             and pin the state. Cannot be reversed without
//                             explicitly clearing is_mutant[i].
// ----------------------------------------------------------------------------
enum class EffectType
{
    INHIBIT,      ///< Multiplicative suppression:  gene *= (1 - strength)
    ACTIVATE,     ///< Multiplicative activation:   gene += strength * (1 - gene)
    KNOCKDOWN,    ///< Pin state low and lock:       gene = clamp(gene - strength, 0,1); lock
    OVEREXPRESS   ///< Pin state high and lock:      gene = clamp(gene + strength, 0,1); lock
};

// ----------------------------------------------------------------------------
// Intervention — a single therapeutic/genetic perturbation on one gene node
//
// The Python EA generates a JSON array of these and writes it to disk.
// The C++ simulation loads it via BooleanNetwork::load_from_json() and
// calls apply_interventions() each timestep after update().
//
// strength is the evolvable parameter in [0, 1]:
//   0.0 → null effect  |  1.0 → maximum effect
// ----------------------------------------------------------------------------
struct Intervention
{
    GeneIndex   gene_index;   ///< Which gene to target (resolved from name by load_from_json)
    EffectType  effect_type;  ///< Type of perturbation
    double      strength;     ///< Effect magnitude in [0.0, 1.0]
    std::string drug_name;    ///< Human-readable label (for logging / EA bookkeeping)
};

// ----------------------------------------------------------------------------
// GeneParams — per-cell parameter overrides loaded from config / EA genome
//
// The EA operates on this struct: it can shift default gene states, modify
// tau values (representing epigenetic or pharmacological changes), or inject
// fixed perturbations. Passed to BooleanNetwork::initialize().
// ----------------------------------------------------------------------------
struct GeneParams
{
    // Additive offsets applied on top of GENE_INFO defaults at initialization.
    // EA can evolve these to represent patient-specific heterogeneity.
    double state_offset[GENE_COUNT]  = {};   ///< Delta applied to initial state [default: 0.0]

    // Multiplicative scale on tau. Values < 1.0 = faster dynamics (e.g., drug
    // that accelerates protein degradation); > 1.0 = slower (stabilization).
    double tau_scale[GENE_COUNT];            ///< Scale factor on time constant [default: 1.0]

    // Drug effect: direct fractional reduction of a gene's target value each step.
    // EA / intervention logic sets these; 0.0 = no drug effect.
    double drug_inhibition[GENE_COUNT] = {}; ///< Fractional inhibition in [0, 1]

    GeneParams();  ///< Constructor: sets all tau_scale to 1.0
};

// ----------------------------------------------------------------------------
// BooleanNetwork — continuous gene state machine for one PhysiCell cell
// ----------------------------------------------------------------------------
class BooleanNetwork
{
public:

    // ---- State ----------------------------------------------------------------

    double    gene_states[GENE_COUNT]; ///< Current continuous gene activity in [0, 1]
    double    tau[GENE_COUNT];         ///< Time constants (minutes) per gene
    bool      is_mutant[GENE_COUNT];   ///< If true, gene is locked and not updated

    /// Structural gene state (GOF/LOF/WT): fixed for the lifetime of a cell.
    /// Set by set_default_genotype() at initialize(); may be overridden by
    /// KNOCKDOWN / OVEREXPRESS interventions. Used by compute_axis_outcome()
    /// to determine dominance voting contributions.
    GeneState genotype[GENE_COUNT];

    // ---- Lifecycle ------------------------------------------------------------

    BooleanNetwork();

    /// Initialize gene states, tau values, and mutant locks for a given cell type.
    /// @param cell_type  TUMOR or STROMA (drives default state selection and mutations)
    /// @param params     EA-evolvable parameter struct (offsets, tau scales, drug effects)
    void initialize(CellType cell_type, const GeneParams& params);

    // ---- PhysiCell synchronization -------------------------------------------

    /// Read gene states from pCell->custom_data[] into gene_states[].
    /// Call this at the start of each phenotype update so the network reflects
    /// any external changes (e.g., drug uptake already applied by PhysiCell).
    void sync_from_cell(PhysiCell::Cell* pCell);

    /// Write gene_states[] back into pCell->custom_data[].
    /// Call this after update() to commit network state to PhysiCell.
    void sync_to_cell(PhysiCell::Cell* pCell) const;

    // ---- Network update ------------------------------------------------------

    /// Run one timestep of the continuous Boolean network.
    ///
    /// Computes target states from Boolean rules conditioned on local
    /// microenvironment signals, then applies the first-order update equation:
    ///   gene[i] += dt * (target[i] - gene[i]) / tau[i]
    ///
    /// Mutant genes are not updated. Drug inhibition in params is applied to
    /// target values before the update step (representing competitive inhibition).
    ///
    /// @param dt          Timestep (minutes); typically PhysiCell dt_phenotype = 6 min
    /// @param oxygen      Local oxygen concentration (mmHg) from microenvironment voxel
    /// @param tgfb_local  Local TGF-beta substrate concentration [dimensionless]
    /// @param shh_local   Local Sonic Hedgehog substrate concentration [dimensionless]
    /// @param drug_local  Local drug substrate concentration [dimensionless]
    /// @param cfg         Runtime thresholds loaded from XML (use ThresholdConfig::load_from_xml())
    void update(double dt,
                double oxygen,
                double tgfb_local,
                double shh_local,
                double drug_local,
                const ThresholdConfig& cfg);

    // ---- Intervention system -------------------------------------------------

    /// Apply a list of interventions to gene_states[] AFTER update().
    ///
    /// Call order each phenotype timestep:
    ///   1. sync_from_cell(pCell)
    ///   2. update(dt, o2, tgfb, shh, drug, cfg)
    ///   3. apply_interventions(interventions)   ← here
    ///   4. sync_to_cell(pCell)
    ///
    /// Reversible effects (INHIBIT / ACTIVATE) are re-applied every call so
    /// their magnitude is maintained while drug is present.
    /// Irreversible effects (KNOCKDOWN / OVEREXPRESS) set is_mutant[i]=true on
    /// first call; subsequent calls are no-ops for those genes because update()
    /// skips mutant genes, but apply_interventions() still pins the state for
    /// safety.
    ///
    /// @param interventions  Vector of Intervention structs (from load_from_json
    ///                       or built programmatically by the EA).
    void apply_interventions(const std::vector<Intervention>& interventions);

    /// Load an intervention list from a JSON file written by the Python EA.
    ///
    /// Expected format:
    /// @code
    /// {
    ///   "interventions": [
    ///     {"gene": "CCND1", "effect": "INHIBIT",    "strength": 0.8, "name": "palbociclib"},
    ///     {"gene": "HAS2",  "effect": "INHIBIT",    "strength": 0.6, "name": "PEGPH20"},
    ///     {"gene": "BCL_XL","effect": "INHIBIT",    "strength": 0.7, "name": "navitoclax"},
    ///     {"gene": "MYC",   "effect": "KNOCKDOWN",  "strength": 0.9, "name": "siRNA_MYC"},
    ///     {"gene": "TP53",  "effect": "OVEREXPRESS","strength": 1.0, "name": "adeno_TP53"}
    ///   ]
    /// }
    /// @endcode
    ///
    /// Gene names are resolved to GeneIndex via GENE_INFO[].name lookup.
    /// Throws std::runtime_error if the file cannot be opened, if JSON is
    /// malformed, if an unknown gene name is encountered, or if effect string
    /// does not match a valid EffectType.
    ///
    /// @param filepath  Absolute or relative path to the JSON file.
    /// @returns         Vector of fully resolved Intervention structs.
    static std::vector<Intervention> load_from_json(const std::string& filepath);

    // ---- Dominance Voting (Gene State Representation & Combination Rules v1.0) ---

    /// Set genotype[] from the canonical default table for this cell type.
    /// Called automatically by initialize(); may be re-called if the EA overrides
    /// a gene's structural state via KNOCKDOWN / OVEREXPRESS.
    void set_default_genotype(CellType cell_type);

    /// Compute the qualitative behavioral axis outcome using dominance voting
    /// over active gene contributions (Part 3–4 of the combination rules document).
    ///
    /// Gene activity thresholds used internally:
    ///   active:  gene_states[i] > 0.50  (gene is functionally ON)
    ///   strong:  gene_states[i] > 0.75  (dominant, e.g., KRAS=1.0)
    ///   low:     gene_states[i] < 0.25  (LOF / silenced)
    ///
    /// Conditional contributions (◆ entries) are evaluated against:
    ///   tgfb_local > TGFB_ACTIVATION_THRESHOLD  (for SMAD4 conditional entries)
    ///   drug_local > 0.01                        (for ABCB1 conditional entry)
    ///
    /// Does NOT apply veto rules — caller must apply veto rules separately
    /// (see tumor_phenotype_update() in tumor_cell.cpp).
    ///
    /// @param axis       Which behavioral axis to evaluate
    /// @param tgfb_local Local TGF-β concentration (for conditional contributions)
    /// @param drug_local Local drug concentration (for ABCB1 conditional)
    /// @returns          Qualitative outcome: NONE → VERY_HIGH
    AxisOutcome compute_axis_outcome(FunctionalAxis axis,
                                     double tgfb_local,
                                     double drug_local) const;

    // ---- Convenience accessors -----------------------------------------------

    /// Return a reference to a gene state by index.
    inline double& operator[](int i)             { return gene_states[i]; }
    inline double  operator[](int i) const        { return gene_states[i]; }
    inline double& operator[](GeneIndex g)        { return gene_states[g]; }
    inline double  operator[](GeneIndex g) const  { return gene_states[g]; }

private:

    // Cached cell type set during initialize(); drives which rule set is applied.
    CellType cell_type_;

    // Cached params set during initialize(); used each update() call.
    GeneParams params_;

    // ---- Rule sets (called by update) ----------------------------------------
    void compute_tumor_targets (double* target, double oxygen,
                                double tgfb_local, double shh_local,
                                double drug_local,
                                const ThresholdConfig& cfg) const;

    void compute_stroma_targets(double* target, double oxygen,
                                double tgfb_local, double shh_local,
                                double drug_local,
                                const ThresholdConfig& cfg) const;

    // ---- Helpers --------------------------------------------------------------

    /// Clamp x to [lo, hi].
    static inline double clamp(double x, double lo, double hi)
    {
        return x < lo ? lo : (x > hi ? hi : x);
    }

    /// Assign default tau values (signaling / TF / structural classification).
    void assign_default_tau();
};

#endif // BOOLEAN_NETWORK_H
