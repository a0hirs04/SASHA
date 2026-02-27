# PROJECT-NORTHSTAR

## Calibration model
- Runtime calibration uses **13 knobs from 8 traits** (`config/tumor_calibration_knobs.json`).
- Active profile defaults to **AsPC-1** and can be switched to **PANC-1** via:
  - XML: `user_parameters/calibration_profile`
  - intervention JSON: `"calibration_profile": "AsPC-1|PANC-1"`

## Partition policy (critical)
- Knobs are partitioned as:
  - `Fixed` (5): tumor-defining, read-only
  - `Observable` (4): measurable state descriptors, read-only
  - `Targetable` (4): EA/drug-touch surface
- Targetable set is strictly:
  - `tgfb_secretion_rate` (1a)
  - `shh_secretion_rate` (1b)
  - `efflux_induction_delay` (7a)
  - `efflux_strength` (7b)
- Partition is enforced in C++ and Python; non-targetable overrides fail early.

## Intervention schema
- Canonical payload:
```json
{
  "calibration_profile": "AsPC-1",
  "knob_interventions": [
    {"knob": "tgfb_secretion_rate", "effect": "INHIBIT", "strength": 0.6, "name": "example"}
  ],
  "drug_delivery": {}
}
```
- One-release legacy bridge is limited to:
  - `TGFB1 -> tgfb_secretion_rate`
  - `SHH -> shh_secretion_rate`
  - `ABCB1 -> efflux_strength`
- Other legacy gene targets are rejected with partition-violation errors.

## EA touch surface
- EA config now uses `targetable_knobs` (not `druggable_genes`).
- Individuals are knob-based with at most one intervention per knob and max 4 interventions.
