# PROJECT-NORTHSTAR

# Stroma World (PDAC Barrier Counterfactual Evolution Engine)

A simulation-first research platform for discovering **stroma-aware** anti-PDAC strategies using a **counterfactual evolution** framework.

> **Current status:** Phase 1 (World Spec) complete and locked. Implementation in progress.

---

## Why this project exists

Pancreatic ductal adenocarcinoma (PDAC) is notoriously resistant to therapy in part because the tumor exists inside a **dense stromal barrier** that limits penetration, alters cell behavior, and creates protected “sanctuaries.” Traditional approaches often optimize therapies in simplified settings and then fail when confronted with barrier dynamics.

This project builds a **minimal-but-faithful PDAC stroma barrier world** and uses **counterfactual evolution** (search over intervention policies under strict safety and realism constraints) to identify strategies that perform well in hostile, PDAC-like microenvironments.

---

## Core idea (plain English)

1. Build a **PDAC-faithful “mini-world”** that reproduces key qualitative PDAC behaviors (barrier formation, drug exclusion, sanctuary formation, and safety tradeoffs).
2. Define a small set of **targetable intervention knobs** representing actions that can plausibly map to real interventions.
3. Use an **evolutionary search** (and later, surrogate acceleration) to find robust strategies that:
   - control tumor burden,
   - avoid catastrophic spread,
   - and respect collateral damage constraints.

This is not “magic AI cures cancer.”  
It is a disciplined engineering pipeline that produces **testable hypotheses and schedules** for validation.

---

## Project roadmap (Phases)

### Phase 0 — Constitution (Mission + Non-negotiables)
Defines hard constraints:
- What counts as “success”
- What is forbidden (no cheating, no omniscience)
- Safety/collateral limits
- Versioning + reproducibility requirements

### Phase 1 — Stroma Barrier World Spec (✅ COMPLETE)
A locked, versioned spec describing:
- Actors (cell types / world objects)
- State variables
- Rules (plain English dynamics)
- Reality anchors + metrics + pass/fail tests
- Module wiring (directed influence graph)
- Assumptions/exclusions and revisit triggers

### Phase 2 — Action Menu Spec (IN PROGRESS / NEXT)
Defines how interventions are represented:
- Action definitions + knobs (intensity, duration, timing, scope)
- Metrics + lab readouts
- Safety costs and legality rules
- Sequencing constraints (PDAC-specific)

### Phase 3 — Counterfactual Evolution Engine
- Strategy/policy search (EA/RL/hybrid)
- Surrogate scoring + truth checks (to make search feasible)
- Anti-cheat stress tests + robustness sweeps

### Phase 4 — Strategy Products
Turn winners into lab-ready playbooks:
- Minimal action sets (2–6 actions)
- Sequencing schedules
- Failure modes + boundary conditions
- Experiment matrices

### Phase 5 — Wet-lab validation
Validate in PDAC-relevant systems (progressive ladder):
- Co-culture → spheroids → organoids → microfluidics → in vivo (later)

### Phase 6 — Calibration loop
Use lab outcomes to refine:
- Parameters
- Missing mechanisms
- Constraint tuning  
(with strict change control + regression tests)

### Phase 7 — Production readiness
Prove it’s a platform:
- Robustness across virtual patients
- Safety audit
- Replayability
- Documentation + versioned release criteria

---

## What’s in this repository

> **NOTE:** repo structure may change as implementation progresses.

Suggested layout:
