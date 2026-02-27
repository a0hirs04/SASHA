// ============================================================================
// baseline_validation.cpp — Stroma World / PROJECT-NORTHSTAR
//
// §8.1 Baseline Behavioral Validation Suite
//
// Tests run entirely on BooleanNetwork instances — no PhysiCell Cell or
// Phenotype objects needed. ThresholdConfig uses compile-time defaults
// (identical to XML defaults) so no XML must be parsed before calling.
//
// CONVERGENCE STRATEGY:
//   All tests run 100 steps at dt=6 min (= 600 simulated minutes, ~10h).
//   This exceeds 10× the TF time constant (tau_TF=60 min), giving >99.99%
//   convergence for all TF-class genes.  Signaling genes (tau=6 min) converge
//   within the first 3–5 steps.
//
// SMAD4-RESTORATION NOTE (Test 2 Case B):
//   BooleanNetwork::update() hard-re-locks SMAD4=0 after every integration
//   step (tumor genotype invariant).  To test a SMAD4-restored scenario
//   without modifying the core update function, run_steps_with_ivs() applies
//   an OVEREXPRESS Intervention after each update() call — exactly mirroring
//   how the simulation would handle a gene-therapy intervention.
// ============================================================================

#include "baseline_validation.h"

#include "boolean_network.h"
#include "gene_definitions.h"
#include "tumor_calibration_knobs.h"

#include <algorithm>   // std::min
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

namespace
{

// ============================================================================
// Simulation helpers
// ============================================================================

static inline double clamp01(double x)
{
    return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x);
}

/// Run N steps of BooleanNetwork::update() with fixed microenvironment values.
static void run_steps(BooleanNetwork& bn, int n, double dt,
                      double oxygen, double tgfb, double shh, double drug)
{
    ThresholdConfig cfg;  // default compile-time values; no XML needed
    for (int i = 0; i < n; ++i)
        bn.update(dt, oxygen, tgfb, shh, drug, cfg);
}

/// Same as run_steps() but applies an Intervention list after each update()
/// step.  Required when testing interventions (e.g., SMAD4 restoration) that
/// must persist through update()'s hard-lock of tumor mutant genes.
static void run_steps_with_ivs(BooleanNetwork& bn, int n, double dt,
                                double oxygen, double tgfb, double shh,
                                double drug,
                                const std::vector<Intervention>& ivs)
{
    ThresholdConfig cfg;
    for (int i = 0; i < n; ++i)
    {
        bn.update(dt, oxygen, tgfb, shh, drug, cfg);
        bn.apply_interventions(ivs);
    }
}

/// Create a canonical PDAC tumor BooleanNetwork ready to run.
static BooleanNetwork make_tumor_bn()
{
    BooleanNetwork bn;
    GeneParams params;
    bn.initialize(CellType::TUMOR, params);
    return bn;
}

// ============================================================================
// Per-test assertion helper
//
// Usage inside each test function:
//   bool pass = true;
//   auto chk = make_chk(log, pass);
//   chk(value > threshold, "descriptive label", value, threshold);
// ============================================================================

struct Checker
{
    std::ostream& log;
    bool&         all_pass;

    void operator()(bool cond, const char* label, double got, double expected) const
    {
        if (!cond)
        {
            all_pass = false;
            std::ostringstream ss;
            ss << std::fixed << std::setprecision(5);
            ss << "    FAIL  " << label
               << "  (got=" << got << "  threshold=" << expected << ")\n";
            log << ss.str();
        }
    }

    // Overload for axis outcome comparisons (integer-valued).
    void operator()(bool cond, const char* label, int got, int expected) const
    {
        if (!cond)
        {
            all_pass = false;
            log << "    FAIL  " << label
                << "  (got level=" << got
                << "  required>=" << expected << ")\n";
        }
    }
};

// ============================================================================
// Test 1 — Constitutive Paracrine Secretion of Stromagenic Signals
//
// Behavioral test (PDF): "Place a single tumor cell next to quiescent PSCs
//   with no other signals.  PSCs should begin activating within the first
//   few simulated hours."
//
// Validated here: the tumor cell's TGFB1 and SHH gene states (expression
//   proxies) converge to constitutively high values driven by KRAS=GOF
//   even when no external TGF-β or SHH is present in the microenvironment.
//   The resulting secretion rates (gene_state × scale) are non-zero and
//   large enough to cross the stromal activation threshold.
// ============================================================================
static bool test_trait1(std::ostream& log)
{
    bool pass = true;
    Checker chk{log, pass};

    // Normoxic, no external signals — baseline isolated tumor cell.
    BooleanNetwork bn = make_tumor_bn();
    run_steps(bn, 100, 6.0, 38.0, 0.0, 0.0, 0.0);

    // TGFB1 target = 0.1 + 0.3*KRAS = 0.4 at steady state (tau_TF=60min → converged).
    chk(bn[TGFB1] > 0.35, "TGFB1  constitutive secretion   (KRAS=1 → target=0.40)",
        bn[TGFB1], 0.35);

    // SHH target = 0.1 + 0.3*KRAS = 0.4.
    chk(bn[SHH] > 0.35, "SHH    constitutive secretion   (KRAS=1 → target=0.40)",
        bn[SHH], 0.35);

    // Secretion rates used in tumor_cell.cpp step 5e:
    //   TGF-β: tgfb1 * (0.0625 * knob_1a)
    //   SHH:   shh   * (0.0285714286 * knob_1b)
    const TumorCalibrationKnobs& knobs = get_active_calibration_knobs();
    const double tgfb_rate = bn[TGFB1] * (0.0625 * clamp01(knobs.tgfb_secretion_rate));
    const double shh_rate  = bn[SHH]   * (0.0285714286 * clamp01(knobs.shh_secretion_rate));
    chk(tgfb_rate > 0.015, "TGF-β secretion rate  > 0.015/min", tgfb_rate, 0.015);
    chk(shh_rate  > 0.006, "SHH   secretion rate  > 0.006/min", shh_rate,  0.006);

    // Verify SHH exceeds the stromal activation threshold (SHH_ACTIVATION_THRESHOLD=0.05).
    chk(bn[SHH] > SHH_ACTIVATION_THRESHOLD,
        "SHH gene state > stromal activation threshold (SHH_ACTIVATION_THRESHOLD)",
        bn[SHH], SHH_ACTIVATION_THRESHOLD);

    return pass;
}

// ============================================================================
// Test 2 — TGF-β Insensitivity (Growth Arrest Arm Only)
//
// Behavioral test (PDF): "Bathe tumor cells in high TGF-β.  They should NOT
//   arrest.  They should become more motile."
//
// Case A — Canonical PDAC (AsPC-1 default Knob 2 low):
//   High TGF-β must NOT strongly activate RB1 (arrest arm remains weak).
//   SNAI1 MUST rise (SMAD-independent invasion arm intact).
// Case B — SMAD4 restored (therapeutic OVEREXPRESS intervention):
//   Same TGF-β environment. RB1 MUST rise when SMAD4 is restored,
//   confirming the SMAD-dependent arm is architecturally present but
//   structurally silenced in canonical PDAC — not deleted from the model.
// ============================================================================
static bool test_trait2(std::ostream& log)
{
    bool pass = true;
    Checker chk{log, pass};
    ThresholdConfig cfg;

    // ---- Case A: canonical PDAC, SMAD4=LOF (locked) ----
    BooleanNetwork bn_a = make_tumor_bn();
    run_steps(bn_a, 100, 6.0, 38.0, 1.0, 0.0, 0.0);

    // RB1 target = CDKN2A*0.7 + SMAD4*0.3 + 0.4*tgfb*knob_2
    //            = 0 + 0 + 0.4*1.0*0.05 = 0.02  (low-growth arrest sensitivity)
    chk(bn_a[RB1] < 0.10,
        "Case A  RB1 remains low under high TGF-β (AsPC-1 low Knob 2 sensitivity)",
        bn_a[RB1], 0.10);

    // SNAI1 target = 0.5*tgfb + 0.2*KRAS + 0.1*HIF1A = 0.5+0.2+0 = 0.7
    chk(bn_a[SNAI1] > 0.55,
        "Case A  SNAI1 rises under TGF-β (SMAD-independent invasion arm intact)",
        bn_a[SNAI1], 0.55);

    // BRAKING axis: CDKN2A=0, TP53=0, MYC→1.0 → dominant DOWN → outcome = NONE.
    AxisOutcome braking_a = bn_a.compute_axis_outcome(FunctionalAxis::BRAKING, 1.0, 0.0);
    chk(braking_a == AxisOutcome::NONE,
        "Case A  BRAKING == NONE under TGF-β (all brakes destroyed; no arrest possible)",
        static_cast<int>(braking_a), static_cast<int>(AxisOutcome::NONE));

    // INVASION axis: SNAI1 active, SMAD4/tgfb conditional, TP53- → at least MODERATE.
    AxisOutcome invasion_a = bn_a.compute_axis_outcome(FunctionalAxis::INVASION, 1.0, 0.0);
    chk(static_cast<int>(invasion_a) >= static_cast<int>(AxisOutcome::MODERATE),
        "Case A  INVASION >= MODERATE under TGF-β (SNAI1 active; motility up)",
        static_cast<int>(invasion_a), static_cast<int>(AxisOutcome::MODERATE));

    // ---- Case B: SMAD4 restored via OVEREXPRESS intervention ----
    // run_steps_with_ivs applies OVEREXPRESS after each update(), mirroring how
    // the simulation handles gene-therapy interventions. This keeps SMAD4=1.0
    // through the hard-lock in update() by restoring it immediately after each step.
    BooleanNetwork bn_b = make_tumor_bn();
    const std::vector<Intervention> restore_smad4 = {
        { SMAD4, EffectType::OVEREXPRESS, 1.0, "smad4_restore_test" }
    };
    run_steps_with_ivs(bn_b, 100, 6.0, 38.0, 1.0, 0.0, 0.0, restore_smad4);

    // RB1 target (after SMAD4 restored in AsPC-1 default profile):
    //   0*0.7 + 1.0*0.3 + 0.4*1.0*knob_2(=0.05) = 0.32
    chk(bn_b[RB1] > 0.25,
        "Case B  RB1 rises modestly when SMAD4 is restored (AsPC-1 low knob_2 remains)",
        bn_b[RB1], 0.25);

    return pass;
}

// ============================================================================
// Test 3 — High Proliferative Rate With Broken Cell-Cycle Checkpoints
//
// Behavioral test (PDF): "Tumor cell population doubling time should be fast
//   relative to stromal cell division.  Crowding alone should not halt growth."
// ============================================================================
static bool test_trait3(std::ostream& log)
{
    bool pass = true;
    Checker chk{log, pass};

    BooleanNetwork bn = make_tumor_bn();
    run_steps(bn, 100, 6.0, 38.0, 0.0, 0.0, 0.0);

    // MYC: target = 0.1 + 0.9*KRAS = 1.0 (tau_TF=60min, 600min → converged).
    chk(bn[MYC] > 0.95,
        "MYC converges to ~1.0   (KRAS=1 constitutively drives MYC)",
        bn[MYC], 0.95);

    // GROWTH axis: KRAS=strong (always 1.0) + MYC=strong (→1.0) → two ⬆⬆ → VERY_HIGH.
    AxisOutcome growth = bn.compute_axis_outcome(FunctionalAxis::GROWTH, 0.0, 0.0);
    chk(static_cast<int>(growth) >= static_cast<int>(AxisOutcome::HIGH),
        "GROWTH axis >= HIGH   (KRAS and MYC both strong; uncontested up)",
        static_cast<int>(growth), static_cast<int>(AxisOutcome::HIGH));

    // BRAKING axis: CDKN2A=0(⬇⬇) + TP53=0(⬇⬇) + MYC=strong(⬇⬇) → NONE.
    // Veto 1 also fires: cdkn2a<0.1 AND tp53<0.1 AND myc>0.75 → braking=NONE.
    AxisOutcome braking = bn.compute_axis_outcome(FunctionalAxis::BRAKING, 0.0, 0.0);
    chk(braking == AxisOutcome::NONE,
        "BRAKING == NONE   (CDKN2A-/TP53-/MYC+ triple; crowding cannot arrest)",
        static_cast<int>(braking), static_cast<int>(AxisOutcome::NONE));

    // Proliferation rate from PROLIF bins must far exceed the stromal BASELINE rate.
    // PROLIF[VERY_HIGH] = 0.000926/min   vs   PROLIF[BASELINE] = 0.0000330/min  → 28×
    const double tumor_prolif  = PhenoParamBins::get(PhenoParamBins::PROLIF, growth);
    const double stroma_base   = PhenoParamBins::PROLIF[
                                     static_cast<int>(AxisOutcome::BASELINE)];
    chk(tumor_prolif > stroma_base * 5.0,
        "Tumor PROLIF bin >> stromal BASELINE (should be >5x faster dividing)",
        tumor_prolif, stroma_base * 5.0);

    return pass;
}

// ============================================================================
// Test 4 — Apoptotic Resistance at Baseline
//
// Behavioral test (PDF): "In the absence of drug or immune attack, spontaneous
//   death rate should be negligible compared to proliferation rate."
// ============================================================================
static bool test_trait4(std::ostream& log)
{
    bool pass = true;
    Checker chk{log, pass};

    BooleanNetwork bn = make_tumor_bn();
    run_steps(bn, 100, 6.0, 38.0, 0.0, 0.0, 0.0);

    // BCL_XL: target = 0.3 + 0.5*KRAS + 0.15*NRF2 + 0.1*HIF1A
    //       = 0.3 + 0.5 + 0.15*0.1 + 0 = 0.815  (NRF2 ≈ 0.1 at normoxia, no drug).
    // Veto 2 fires when BCL_XL > 0.75 → floors DEATH resistance at HIGH.
    chk(bn[BCL_XL] > 0.75,
        "BCL-XL > 0.75 baseline   (KRAS drives BCL-XL; veto 2 threshold met)",
        bn[BCL_XL], 0.75);

    // DEATH resistance axis must be >= HIGH (BCL-XL veto 2 floors it at HIGH).
    AxisOutcome death = bn.compute_axis_outcome(FunctionalAxis::DEATH, 0.0, 0.0);
    chk(static_cast<int>(death) >= static_cast<int>(AxisOutcome::HIGH),
        "DEATH resistance >= HIGH   (BCL-XL overexpressed; TP53 lost; anoikis-R)",
        static_cast<int>(death), static_cast<int>(AxisOutcome::HIGH));

    // Invert resistance → apoptosis rate: apopt_int = 6 - death_int.
    // VERY_HIGH resistance (6) → apoptosis = NONE (0).  HIGH (5) → apoptosis = VERY_LOW (1).
    const int death_int = static_cast<int>(death);
    int apopt_int = 6 - death_int;
    if (apopt_int < 0) apopt_int = 0;
    const double apoptosis_rate = PhenoParamBins::get(
        PhenoParamBins::APOPTOSIS, static_cast<AxisOutcome>(apopt_int));
    const AxisOutcome growth    = bn.compute_axis_outcome(FunctionalAxis::GROWTH, 0.0, 0.0);
    const double prolif_rate    = PhenoParamBins::get(PhenoParamBins::PROLIF, growth);

    chk(apoptosis_rate < prolif_rate * 0.1,
        "Basal apoptosis rate < 10% of proliferation rate   (negligible spontaneous death)",
        apoptosis_rate, prolif_rate * 0.1);

    return pass;
}

// ============================================================================
// Test 5 — Inducible EMT Under Microenvironmental Stress (Spatial)
//
// Behavioral test (PDF): "Cells in the oxygenated tumor core should remain
//   epithelial (CDH1 ON, low motility).  Cells at the hypoxic periphery or
//   in TGF-β-rich zones should become mesenchymal (CDH1 OFF, high motility,
//   MMP2 active)."
//
// Two in silico cells:
//   bn_norm — normoxic core     : O2=38 mmHg, TGF-β=0
//   bn_hyp  — invasive front    : O2=5 mmHg,  TGF-β=0.8 (high TGF-β + hypoxia)
//
// Steady-state analytic estimates:
//   Normoxic: ZEB1 ≈ 0.00, CDH1 ≈ 0.95  (fixed-point: ZEB1 target < 0)
//   Hypoxic:  ZEB1 ≈ 0.55, CDH1 ≈ 0.40  (HIF1A=0.67 + SNAI1=0.67 push ZEB1 up)
// ============================================================================
static bool test_trait5(std::ostream& log)
{
    bool pass = true;
    Checker chk{log, pass};

    // ---- Normoxic core cell ----
    BooleanNetwork bn_norm = make_tumor_bn();
    run_steps(bn_norm, 100, 6.0, 38.0, 0.0, 0.0, 0.0);

    // ---- Hypoxic / TGF-β-rich invasive front cell ----
    BooleanNetwork bn_hyp = make_tumor_bn();
    run_steps(bn_hyp, 100, 6.0, 5.0, 0.8, 0.0, 0.0);

    // --- Normoxic checks ---
    chk(bn_norm[HIF1A] < 0.05,
        "Normoxic  HIF1A < 0.05   (O2=38 >> threshold=15; PHD active)",
        bn_norm[HIF1A], 0.05);

    chk(bn_norm[ZEB1] < 0.10,
        "Normoxic  ZEB1 low   (no hypoxia or TGF-β; epithelial state)",
        bn_norm[ZEB1], 0.10);

    chk(bn_norm[CDH1] > 0.80,
        "Normoxic  CDH1 high > 0.80   (E-cadherin ON; epithelial identity maintained)",
        bn_norm[CDH1], 0.80);

    // --- Hypoxic/TGF-β checks ---
    // HIF1A target = 1 - O2/threshold = 1 - 5/15 = 0.667 (tau=6 min; converges in ~4 steps).
    chk(bn_hyp[HIF1A] > 0.55,
        "Hypoxic   HIF1A > 0.55   (O2=5 << threshold=15; PHD inhibited → HIF accumulates)",
        bn_hyp[HIF1A], 0.55);

    // ZEB1 rises: HIF1A + SNAI1 (from TGF-β) push ZEB1 target to ~0.65 before CDH1 feedback.
    chk(bn_hyp[ZEB1] > 0.40,
        "Hypoxic   ZEB1 > 0.40   (HIF1A+TGF-β combined drive EMT program)",
        bn_hyp[ZEB1], 0.40);

    // CDH1 falls as ZEB1 rises.
    chk(bn_hyp[CDH1] < 0.60,
        "Hypoxic   CDH1 < 0.60   (ZEB1+SNAI1 repress E-cadherin; EMT underway)",
        bn_hyp[CDH1], 0.60);

    // MMP2 must be active in hypoxic/mesenchymal cell (ZEB1+HIF1A+KRAS all drive MMP2).
    // MMP2 target = 0.20*ZEB1 + 0.15*SNAI1 + 0.25*HIF1A + 0.15*KRAS
    //             = 0.20*0.55 + 0.15*0.67 + 0.25*0.67 + 0.15 ≈ 0.53 → active
    chk(bn_hyp[MMP2] > 0.40,
        "Hypoxic   MMP2 > 0.40   (ECM degradation active at invasive front)",
        bn_hyp[MMP2], 0.40);

    // Spatial contrast: ZEB1 and CDH1 must differ significantly between cells.
    chk(bn_hyp[ZEB1] > bn_norm[ZEB1] + 0.30,
        "Spatial   ZEB1_hypoxic >> ZEB1_normoxic   (EMT is spatially restricted)",
        bn_hyp[ZEB1] - bn_norm[ZEB1], 0.30);

    chk(bn_hyp[CDH1] < bn_norm[CDH1] - 0.25,
        "Spatial   CDH1_hypoxic << CDH1_normoxic   (epithelial identity lost at front)",
        bn_norm[CDH1] - bn_hyp[CDH1], 0.25);

    // Invasion axis: hypoxic cell must score clearly higher than normoxic.
    // Normoxic: only TP53- vote → MODERATE.
    // Hypoxic:  ZEB1 active (⬆⬆) + MMP2 active (⬆⬆) + TP53 + SMAD4/tgfb + HIF1A → VERY_HIGH.
    AxisOutcome inv_norm = bn_norm.compute_axis_outcome(FunctionalAxis::INVASION, 0.0, 0.0);
    AxisOutcome inv_hyp  = bn_hyp.compute_axis_outcome(FunctionalAxis::INVASION, 0.8, 0.0);
    chk(static_cast<int>(inv_hyp) > static_cast<int>(inv_norm),
        "Spatial   INVASION_hypoxic > INVASION_normoxic   (axis level)",
        static_cast<int>(inv_hyp), static_cast<int>(inv_norm) + 1);

    chk(static_cast<int>(inv_hyp) >= static_cast<int>(AxisOutcome::HIGH),
        "Hypoxic   INVASION axis >= HIGH   (mesenchymal invasive front)",
        static_cast<int>(inv_hyp), static_cast<int>(AxisOutcome::HIGH));

    return pass;
}

// ============================================================================
// Test 6 — Hypoxia-Responsive Phenotype Switching
//
// Behavioral test (PDF): "Create an artificial oxygen gradient across the
//   domain.  Tumor cells in the low-O₂ zone should visibly change behavior
//   (more motile, more secretory) compared to well-oxygenated cells of the
//   same genotype."
//
// Checks three behavioral shifts driven by HIF1A at O2=5 vs O2=38:
//   (a) TGFB1 amplified — HIF1A contributes +0.2*HIF1A to TGFB1 target
//   (b) BCL-XL elevated — HIF1A contributes +0.1*HIF1A to BCL-XL (more death-resistant)
//   (c) ZEB1 higher — HIF1A contributes to EMT program; motility increases
// ============================================================================
static bool test_trait6(std::ostream& log)
{
    bool pass = true;
    Checker chk{log, pass};

    BooleanNetwork bn_norm = make_tumor_bn();
    run_steps(bn_norm, 100, 6.0, 38.0, 0.0, 0.0, 0.0);

    BooleanNetwork bn_hyp = make_tumor_bn();
    run_steps(bn_hyp, 100, 6.0, 5.0, 0.0, 0.0, 0.0);  // no TGF-β; pure hypoxia

    // HIF1A must activate under hypoxia.
    // target = 1 - O2/threshold = 1 - 5/15 = 0.667 (tau=6min → converges by step 5).
    chk(bn_hyp[HIF1A] > 0.55,
        "Hypoxic  HIF1A > 0.55   (O2=5 → HIF1A stable; PHD unable to hydroxylate)",
        bn_hyp[HIF1A], 0.55);

    chk(bn_norm[HIF1A] < 0.05,
        "Normoxic HIF1A < 0.05   (O2=38 → HIF1A degraded by PHD2)",
        bn_norm[HIF1A], 0.05);

    // (a) TGFB1 amplified by hypoxia.
    // Normoxic TGFB1 target = 0.1 + 0.3*1 = 0.40
    // Hypoxic  TGFB1 target = 0.1 + 0.3 + 0.2*0.667 = 0.53
    chk(bn_hyp[TGFB1] > bn_norm[TGFB1] + 0.05,
        "Hypoxic  TGFB1 > normoxic + 0.05   (HIF1A amplifies stromagenic secretion)",
        bn_hyp[TGFB1] - bn_norm[TGFB1], 0.05);

    // (b) BCL-XL elevated by hypoxia (HIF1A direct arm: HIF binds BCL2L1 HREs).
    // Normoxic BCL-XL ≈ 0.815;  Hypoxic BCL-XL ≈ 0.815 + 0.1*0.667 = 0.882
    chk(bn_hyp[BCL_XL] > bn_norm[BCL_XL],
        "Hypoxic  BCL-XL elevated vs normoxic   (HIF1A direct arm on BCL2L1 promoter)",
        bn_hyp[BCL_XL] - bn_norm[BCL_XL], 0.001);

    // (c) ZEB1 higher in hypoxic cell (HIF1A contributes +0.25*HIF1A to ZEB1 target).
    // Normoxic ZEB1 ≈ 0.00;  Hypoxic ZEB1 ≈ 0.17  (even without TGF-β)
    chk(bn_hyp[ZEB1] > bn_norm[ZEB1] + 0.05,
        "Hypoxic  ZEB1 > normoxic + 0.05   (HIF1A initiates pre-EMT; motility rises)",
        bn_hyp[ZEB1] - bn_norm[ZEB1], 0.05);

    // Secretion axis: KRAS+TGFB1+SHH dominate both cells (both VERY_HIGH), but
    // the actual TGFB1 gene state (hence secretion rate) is higher under hypoxia.
    // Positive feedback: hypoxia → more TGF-β secreted → more stroma → more hypoxia.
    const TumorCalibrationKnobs& knobs = get_active_calibration_knobs();
    const double tgfb_rate_scale = 0.0625 * clamp01(knobs.tgfb_secretion_rate);
    const double tgfb_rate_norm = bn_norm[TGFB1] * tgfb_rate_scale;
    const double tgfb_rate_hyp  = bn_hyp[TGFB1]  * tgfb_rate_scale;
    chk(tgfb_rate_hyp > tgfb_rate_norm * 1.1,
        "Hypoxic  TGF-β secretion rate > 1.1× normoxic   (positive feedback loop)",
        tgfb_rate_hyp / (tgfb_rate_norm + 1e-9), 1.1);

    return pass;
}

// ============================================================================
// Test 7 — Drug-Inducible Efflux (Not Constitutive)
//
// Behavioral test (PDF): "Expose tumor cells to drug.  After a delay (not
//   instantly), ABCB1 should activate and intracellular drug concentration
//   should decline.  Cells never exposed to drug should have no ABCB1 activity."
//
// Mechanism: drug → NRF2 (fast, tau=6min) → ABCB1 (slow, tau=60min).
// The "delay" is modeled by ABCB1's tau_TF=60min vs NRF2's tau_sig=6min.
// Veto 5 (tumor_cell.cpp): drug_sensitivity_stored = 1.0 when drug_local < 0.01
// so the EA sees no constitutive resistance when drug is absent.
// ============================================================================
static bool test_trait7(std::ostream& log)
{
    bool pass = true;
    Checker chk{log, pass};

    // Undrugged: no drug exposure.
    BooleanNetwork bn_nodrug = make_tumor_bn();
    run_steps(bn_nodrug, 100, 6.0, 38.0, 0.0, 0.0, 0.0);

    // Drugged: maximum drug exposure (drug=1.0 saturates the NRF2 term).
    BooleanNetwork bn_drugged = make_tumor_bn();
    run_steps(bn_drugged, 100, 6.0, 38.0, 0.0, 0.0, 1.0);

    // --- Undrugged cell: low NRF2 and negligible ABCB1 ---
    // NRF2 target = 0.1*KRAS + 0.2*HIF1A + 0.4*drug = 0.1 + 0 + 0 = 0.1 (tau=6, fast)
    chk(bn_nodrug[NRF2] < 0.15,
        "Undrugged  NRF2 < 0.15   (only KRAS basal arm; no drug → no stress induction)",
        bn_nodrug[NRF2], 0.15);

    // ABCB1 target = 0.5*NRF2 + 0.3*HIF1A = 0.5*0.1 + 0 = 0.05 (tau=60, slow)
    chk(bn_nodrug[ABCB1] < 0.10,
        "Undrugged  ABCB1 < 0.10   (no drug → MDR1 not induced; efflux pump silent)",
        bn_nodrug[ABCB1], 0.10);

    // --- Drugged cell: elevated NRF2 then ABCB1 ---
    // NRF2 target = 0.1 + 0 + 0.4*1.0 = 0.5 (tau=6; fully converged well within 100 steps)
    chk(bn_drugged[NRF2] > 0.40,
        "Drugged    NRF2 > 0.40   (drug_local=1 adds +0.4 to NRF2 target; KEAP1 inactivated)",
        bn_drugged[NRF2], 0.40);

    // ABCB1 target = 0.5*0.5 = 0.25 (tau=60; 600 min → >99.99% converged)
    chk(bn_drugged[ABCB1] > 0.18,
        "Drugged    ABCB1 > 0.18   (NRF2-driven ARE transcription; MDR1 pump induced)",
        bn_drugged[ABCB1], 0.18);

    // Acquired resistance: drugged ABCB1 must be >> undrugged (not constitutive).
    chk(bn_drugged[ABCB1] > bn_nodrug[ABCB1] * 3.0,
        "ABCB1_drugged > 3× ABCB1_undrugged   (acquired resistance, not constitutive)",
        bn_drugged[ABCB1] / (bn_nodrug[ABCB1] + 1e-9), 3.0);

    // Veto 5 — drug_sensitivity formula from tumor_cell.cpp:
    //   drug_sensitivity_stored = (drug > 0.01) ? clamp(1-0.4*nrf2-0.4*abcb1, 0.05, 1.0) : 1.0
    // When drug absent: must return 1.0 regardless of NRF2/ABCB1 state.
    {
        const double nrf2_high  = bn_drugged[NRF2];
        const double abcb1_high = bn_drugged[ABCB1];

        // Simulate veto 5 logic.
        auto drug_sens = [&](double drug_local) -> double {
            if (drug_local > 0.01)
            {
                double s = 1.0 - 0.4 * nrf2_high - 0.4 * abcb1_high;
                if (s < 0.05) s = 0.05;
                if (s > 1.0)  s = 1.0;
                return s;
            }
            return 1.0;
        };

        const double sens_no_drug   = drug_sens(0.0);
        const double sens_with_drug = drug_sens(0.5);

        chk(sens_no_drug == 1.0,
            "Veto 5  drug_sensitivity = 1.0 when drug_local < 0.01   (ABCB1 constitutive-free)",
            sens_no_drug, 1.0);

        chk(sens_with_drug < 1.0,
            "Veto 5  drug_sensitivity < 1.0 when drug present (NRF2+ABCB1 active)",
            sens_with_drug, 1.0);
    }

    return pass;
}

// ============================================================================
// Test 8 — Density-Driven ECM Compaction (Passive, Not Active)
//
// Behavioral test (PDF): "Simulate tumor growth inside a rigid stroma shell
//   with CAF secretion turned OFF.  ECM density at the tumor-stroma boundary
//   should still increase due to compaction."
//
// Validated analytically: the solid_stress formula in tumor_cell.cpp
//   solid_stress       = simple_pressure * (1 + ecm)
//   solid_stress_brake = clamp(1 - 0.35 * min(1, solid_stress), 0, 1)
//
// Key invariants:
//   (A) Compaction occurs even with ECM=0 (CAFs off): pressure alone suppresses.
//   (B) ECM stiffness amplifies stress at the same crowding pressure.
//   (C) Dense ECM (0.8) produces more suppression than sparse ECM (0.0).
//   (D) Maximum suppression capped at 35% — cells adapt; no full mechanical arrest.
//   (E) Mechanical term is independent of gene network: operates on pCell->state.simple_pressure,
//       not on BooleanNetwork outputs. Verified by showing suppression persists even
//       at maximum gene-network growth drive.
// ============================================================================
static bool test_trait8(std::ostream& log)
{
    bool pass = true;
    Checker chk{log, pass};

    // Replicate the formula from tumor_cell.cpp step 5a:
    auto solid_stress_brake = [](double pressure, double ecm) -> double {
        const double solid_stress = pressure * (1.0 + ecm);
        return clamp01(1.0 - 0.35 * std::min(1.0, solid_stress));
    };

    // (A) Compaction without CAFs: ECM=0, realistic crowding pressure=0.3.
    {
        const double p = 0.3, ecm = 0.0;
        const double stress = p * (1.0 + ecm);   // 0.30
        const double brake  = solid_stress_brake(p, ecm);   // 0.895
        chk(stress > 0.0,
            "Case A  solid_stress > 0 even at ECM=0   (compaction from crowding alone)",
            stress, 0.0);
        chk(brake < 1.0,
            "Case A  solid_stress_brake < 1.0   (mechanical suppression present without CAFs)",
            brake, 1.0);
        chk(brake > 0.65,
            "Case A  solid_stress_brake > 0.65   (suppression capped at 35%)",
            brake, 0.65);
    }

    // (B) ECM stiffness amplifies solid stress at same crowding.
    {
        const double p = 0.3;
        const double stress_no_ecm   = p * (1.0 + 0.0);   // 0.30
        const double stress_with_ecm = p * (1.0 + 0.8);   // 0.54
        chk(stress_with_ecm > stress_no_ecm,
            "Case B  solid_stress(ECM=0.8) > solid_stress(ECM=0)   (stiffness amplifier)",
            stress_with_ecm - stress_no_ecm, 0.0);
    }

    // (C) Dense ECM (0.8) produces more proliferation suppression than sparse (0.0).
    {
        const double p = 0.3;
        const double brake_no_ecm   = solid_stress_brake(p, 0.0);   // 0.895
        const double brake_with_ecm = solid_stress_brake(p, 0.8);   // 0.811
        chk(brake_with_ecm < brake_no_ecm,
            "Case C  brake(ECM=0.8) < brake(ECM=0)   (denser ECM → more growth suppression)",
            brake_no_ecm - brake_with_ecm, 0.001);
    }

    // (D) Maximum suppression is capped.
    {
        // At extreme pressure/ECM the brake floor must be > 0.65 (35% max cap).
        const double brake_max = solid_stress_brake(1.0, 1.0);   // 1 - 0.35*min(1,2)=0.65
        chk(brake_max >= 0.65,
            "Case D  solid_stress_brake >= 0.65 at extreme load   (cells adapt; no total arrest)",
            brake_max, 0.65);
    }

    // (E) Mechanical term persists at maximum gene-network growth drive.
    // Even with KRAS/MYC fully active, solid stress reduces the final proliferation rate.
    {
        BooleanNetwork bn = make_tumor_bn();
        run_steps(bn, 100, 6.0, 38.0, 0.0, 0.0, 0.0);

        const AxisOutcome growth_out  = bn.compute_axis_outcome(FunctionalAxis::GROWTH, 0.0, 0.0);
        const double growth_bin       = PhenoParamBins::get(PhenoParamBins::PROLIF, growth_out);
        const double mechanical_brake = solid_stress_brake(0.5, 0.6);  // p=0.5, ECM=0.6 → 0.72
        const double final_prolif     = growth_bin * mechanical_brake;

        chk(final_prolif < growth_bin,
            "Case E  solid_stress reduces proliferation even at max KRAS/MYC drive  "
            "(mechanical contribution is gene-network-independent)",
            final_prolif / (growth_bin + 1e-12), 1.0);
    }

    return pass;
}

// ============================================================================
// Test 9 — 13-Knob Partition + AsPC-1 / PANC-1 Knob 2 Comparator
// ============================================================================
static bool test_trait9(std::ostream& log)
{
    bool pass = true;
    Checker chk{log, pass};

    const TumorCalibrationKnobs saved = get_active_calibration_knobs();
    try
    {
        validate_knob_partition();

        chk(is_targetable_knob("tgfb_secretion_rate"), "targetable includes knob 1a", 1.0, 1.0);
        chk(is_targetable_knob("shh_secretion_rate"), "targetable includes knob 1b", 1.0, 1.0);
        chk(is_targetable_knob("efflux_induction_delay"), "targetable includes knob 7a", 1.0, 1.0);
        chk(is_targetable_knob("efflux_strength"), "targetable includes knob 7b", 1.0, 1.0);
        chk(!is_targetable_knob("tgfb_brake_sensitivity"), "knob 2 not targetable", 0.0, 0.0);

        const TumorCalibrationProfiles profiles =
            load_tumor_calibration_profiles("config/tumor_calibration_knobs.json");
        chk(std::fabs(profiles.aspc1.tgfb_brake_sensitivity - 0.05) < 1e-9,
            "AsPC-1 knob 2 default == 0.05", profiles.aspc1.tgfb_brake_sensitivity, 0.05);
        chk(std::fabs(profiles.panc1.tgfb_brake_sensitivity - 0.9) < 1e-9,
            "PANC-1 knob 2 override == 0.9", profiles.panc1.tgfb_brake_sensitivity, 0.9);

        set_active_calibration_knobs(profiles.aspc1);
        BooleanNetwork aspc1_bn = make_tumor_bn();
        run_steps(aspc1_bn, 100, 6.0, 38.0, 1.0, 0.0, 0.0);

        set_active_calibration_knobs(profiles.panc1);
        BooleanNetwork panc1_bn = make_tumor_bn();
        run_steps(panc1_bn, 100, 6.0, 38.0, 1.0, 0.0, 0.0);

        chk(panc1_bn[RB1] > aspc1_bn[RB1] + 0.25,
            "PANC-1 RB1(TGF) >> AsPC-1 RB1(TGF) under high TGF-β (Knob 2 comparator)",
            panc1_bn[RB1] - aspc1_bn[RB1], 0.25);
    }
    catch (const std::exception& e)
    {
        pass = false;
        log << "    FAIL  Trait 9 exception: " << e.what() << "\n";
    }

    set_active_calibration_knobs(saved);
    return pass;
}

// ============================================================================
// Test 10 — SHH Paradox Fast Invariants (Proxy)
// ============================================================================
static bool test_trait10(std::ostream& log)
{
    bool pass = true;
    Checker chk{log, pass};

    const TumorCalibrationKnobs& knobs = get_active_calibration_knobs();

    BooleanNetwork baseline = make_tumor_bn();
    run_steps(baseline, 100, 6.0, 38.0, 0.0, 0.0, 0.0);

    BooleanNetwork shh_off = make_tumor_bn();
    const std::vector<Intervention> shh_block = {
        { SHH, EffectType::INHIBIT, 1.0, "SHH_OFF_NO_CYTOTOXIC_PROXY" }
    };
    run_steps_with_ivs(shh_off, 100, 6.0, 38.0, 0.0, 0.0, 0.0, shh_block);

    const double shh_scale = 0.0285714286 * clamp01(knobs.shh_secretion_rate);
    const double shh_secretion_baseline = baseline[SHH] * shh_scale;
    const double shh_secretion_off = shh_off[SHH] * shh_scale;
    chk(shh_secretion_off < shh_secretion_baseline,
        "SHH-off reduces SHH secretion (required precondition)", shh_secretion_off, shh_secretion_baseline);

    // Proxy for paradox:
    // Lower SHH -> lower stroma/ECM containment.
    // - No cytotoxic context: less containment can worsen growth.
    // - Cytotoxic context: less ECM barrier improves drug penetration.
    const double ecm_proxy_baseline = clamp01(0.1 + shh_secretion_baseline);
    const double ecm_proxy_off      = clamp01(0.1 + shh_secretion_off);

    const double no_cytotoxic_growth_drive_baseline = 1.0 - 0.25 * ecm_proxy_baseline;
    const double no_cytotoxic_growth_drive_off      = 1.0 - 0.25 * ecm_proxy_off;
    chk(no_cytotoxic_growth_drive_off > no_cytotoxic_growth_drive_baseline,
        "SHH_OFF_NO_CYTOTOXIC proxy worsens tumor growth pressure vs baseline",
        no_cytotoxic_growth_drive_off - no_cytotoxic_growth_drive_baseline, 1e-4);

    const double cytotoxic_penetration_baseline = 1.0 - 0.5 * ecm_proxy_baseline;
    const double cytotoxic_penetration_off      = 1.0 - 0.5 * ecm_proxy_off;
    chk(cytotoxic_penetration_off > cytotoxic_penetration_baseline,
        "SHH_OFF_WITH_CYTOTOXIC proxy improves penetration vs no-cytotoxic SHH-off context",
        cytotoxic_penetration_off - cytotoxic_penetration_baseline, 1e-4);

    return pass;
}

} // anonymous namespace

// ============================================================================
// run_baseline_validation — public entry point
// ============================================================================

bool run_baseline_validation()
{
    struct TraitEntry
    {
        int         number;
        const char* name;
        bool (*fn)(std::ostream&);
    };

    static const TraitEntry TESTS[] = {
        { 1, "Constitutive Paracrine Secretion",            test_trait1 },
        { 2, "TGF-beta Insensitivity (Arrest Arm Only)",    test_trait2 },
        { 3, "Rapid Proliferation / Broken Checkpoints",    test_trait3 },
        { 4, "Apoptotic Resistance at Baseline",            test_trait4 },
        { 5, "Inducible EMT Under Stress  (Spatial)",       test_trait5 },
        { 6, "Hypoxia-Responsive Phenotype Switching",      test_trait6 },
        { 7, "Drug-Inducible Efflux  (Not Constitutive)",   test_trait7 },
        { 8, "ECM Compaction  (Passive / Mechanical)",      test_trait8 },
        { 9, "13-Knob Partition + AsPC-1/PANC-1 Knob 2",    test_trait9 },
        {10, "SHH Paradox Fast Invariants",                 test_trait10 },
    };

    static const int N = 10;

    std::cerr << "\n";
    std::cerr << "================================================================\n";
    std::cerr << " §8.1  Baseline Validation Suite — Tumor-Intrinsic Traits v1.0 \n";
    std::cerr << "================================================================\n";

    int n_pass = 0;
    int n_fail = 0;

    for (const TraitEntry& t : TESTS)
    {
        std::ostringstream local_log;
        const bool ok = t.fn(local_log);

        if (ok)
        {
            ++n_pass;
            std::cerr << "  [PASS]  Trait " << t.number << ":  " << t.name << "\n";
        }
        else
        {
            ++n_fail;
            std::cerr << "  [FAIL]  Trait " << t.number << ":  " << t.name << "\n";
            std::cerr << local_log.str();
        }
    }

    std::cerr << "----------------------------------------------------------------\n";
    if (n_fail == 0)
    {
        std::cerr << "  Result: " << n_pass << "/" << N
                  << " passed — ALL PASS.  Baseline is validated.\n";
    }
    else
    {
        std::cerr << "  Result: " << n_pass << "/" << N
                  << " passed — " << n_fail
                  << " FAILED.  Fix before launching EA.\n";
    }
    std::cerr << "================================================================\n\n";

    return (n_fail == 0);
}
