#ifndef BASELINE_VALIDATION_H
#define BASELINE_VALIDATION_H

// ============================================================================
// baseline_validation.h — Stroma World / PROJECT-NORTHSTAR
//
// §8.1 Baseline Behavioral Validation Suite
//       "Tumor-Intrinsic Baseline Traits v1.0"
//
// PURPOSE:
//   Verifies that the BooleanNetwork gene-network logic faithfully implements
//   all eight tumor-intrinsic baseline traits before any EA evolution begins.
//   Also enforces 13-knob partition invariants (Fixed/Observable/Targetable)
//   and selected fast calibration checks.
//
//   Each trait maps directly to one or more BooleanNetwork behaviors that
//   can be checked analytically (no full PhysiCell simulation required).
//   Tests operate on BooleanNetwork instances with default ThresholdConfig
//   values — no XML parsing, no Cell/Phenotype objects needed.
//
// CALL SITE:
//   // In setup_tissue() immediately after cell types are registered:
//   if (!run_baseline_validation()) {
//       std::cerr << "[FATAL] Baseline traits failed — do not proceed.\n";
//       exit(EXIT_FAILURE);
//   }
//
// TRAITS / INVARIANTS TESTED:
//   1. Constitutive Paracrine Secretion (TGFB1, SHH driven by KRAS=GOF)
//   2. TGF-β Insensitivity (low Knob 2 growth-arrest sensitivity; invasion arm intact)
//   3. Rapid Proliferation / Broken Checkpoints (KRAS+/MYC+, braking=NONE)
//   4. Apoptotic Resistance at Baseline (BCL-XL+/TP53-; apoptosis << proliferation)
//   5. Inducible EMT (hypoxic/TGF-β periphery mesenchymal; normoxic core epithelial)
//   6. Hypoxia-Responsive Phenotype Switching (HIF1A activates; secretion amplified)
//   7. Drug-Inducible Efflux — not constitutive (NRF2→ABCB1 only when drug present)
//   8. ECM Compaction — mechanical not signaling (solid_stress formula verified)
//   9. 13-knob partition constraints and AsPC-1 / PANC-1 Knob 2 comparator
//  10. SHH paradox fast invariants (SHH-off no-cytotoxic vs SHH-off with cytotoxic)
// ============================================================================

/// Run all 8 baseline behavioral tests.
///
/// Writes a PASS/FAIL summary table to std::cerr.
/// @returns  true if all tests pass; false if any test fails.
bool run_baseline_validation();

#endif // BASELINE_VALIDATION_H
