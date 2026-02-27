#ifndef GENE_DEFINITIONS_H
#define GENE_DEFINITIONS_H

// ============================================================================
// gene_definitions.h — Stroma World / PROJECT-NORTHSTAR
//
// Defines all gene-level metadata for the 2D PDAC tumor-stroma simulation:
//   - GeneIndex enum (21 genes + GENE_COUNT sentinel)
//   - GeneInfo struct (name, cell type, functional axis, druggability, defaults)
//   - GENE_INFO[] static metadata array
//   - Phenotype mapping threshold constants
//
// Usage:
//   #include "gene_definitions.h"
//   double val = pCell->custom_data[KRAS];       // index lookup
//   auto&  info = GENE_INFO[ZEB1];               // metadata lookup
//   if (o2 < HYPOXIA_O2_THRESHOLD) { ... }
// ============================================================================

#include <array>
#include <string_view>

// ----------------------------------------------------------------------------
// 1. GENE INDEX ENUM
//    Provides named integer indices into PhysiCell custom_data arrays.
//    Ordering must match the sequence declared in PhysiCell_settings.xml
//    custom_data blocks for tumor_cell and stromal_cell.
//    GENE_COUNT is the sentinel used to size static arrays.
// ----------------------------------------------------------------------------
enum GeneIndex : int
{
    // ---- Tumor oncogenes & signaling ----
    KRAS   = 0,   ///< KRAS G12D/V; constitutively activates RAS/MAPK, PI3K
    MYC    = 1,   ///< MYC transcription factor; drives proliferation, metabolism
    CCND1  = 2,   ///< Cyclin D1; G1/S accelerator; pairs with CDK4 freed by CDKN2A LOF

    // ---- Tumor suppressors / pro-death ----
    TP53   = 3,   ///< p53 guardian of genome; lost in ~75% PDAC
    BCL_XL = 4,   ///< Anti-apoptotic BCL-2 family; navitoclax target
    SNAI1  = 5,   ///< Snail EMT initiator TF; fast CDH1 repressor; first responder to TGF-b

    // ---- Cell-cycle brakes ----
    CDKN2A = 6,   ///< p16/INK4A CDK4/6 inhibitor; deleted in ~95% PDAC
    SMAD4  = 7,   ///< SMAD4 TGF-beta canonical effector; lost in ~55% PDAC
    RB1    = 8,   ///< Retinoblastoma protein; G1/S checkpoint

    // ---- Invasion / EMT ----
    ZEB1   = 9,   ///< ZEB1 EMT transcription factor; represses CDH1
    CDH1   = 10,  ///< E-cadherin; epithelial adhesion; lost during EMT
    MMP2   = 11,  ///< Matrix metalloproteinase-2; ECM degradation (tumor+stroma)

    // ---- Stroma activation markers ----
    ACTA2  = 12,  ///< Alpha-SMA; canonical CAF activation marker

    // ---- Secreted factors (gene expression proxies) ----
    TGFB1  = 13,  ///< TGF-beta1 expression; drives stroma activation (tumor+stroma)
    SHH    = 14,  ///< Sonic Hedgehog ligand expression; tumor->stroma signal
    GLI1   = 15,  ///< GLI1 Hedgehog TF; activated in stroma by SHH
    HAS2   = 16,  ///< Hyaluronan synthase 2; ECM production; PEGPH20 target
    COL1A1 = 17,  ///< Collagen I alpha-1; fibrous ECM deposition

    // ---- Stress response ----
    HIF1A  = 18,  ///< HIF-1alpha; hypoxia master regulator (tumor+stroma)
    NRF2   = 19,  ///< NRF2; oxidative stress / antioxidant response
    ABCB1  = 20,  ///< MDR1/P-gp drug efflux pump; acquired drug resistance

    GENE_COUNT = 21  ///< Sentinel — always last; used to size GENE_INFO[]
};

// ----------------------------------------------------------------------------
// 2. SUPPORTING ENUMS
// ----------------------------------------------------------------------------

/// Which cell population expresses / uses this gene
enum class CellType
{
    TUMOR,   ///< Expressed only in tumor cells
    STROMA,  ///< Expressed only in stromal cells (CAFs/PSCs)
    BOTH     ///< Expressed in both cell types
};

/// Broad functional role in the PhysiCell phenotype model
enum class FunctionalAxis
{
    GROWTH,    ///< Drives proliferation (cycle rate upregulation)
    DEATH,     ///< Controls apoptosis/necrosis sensitivity
    BRAKING,   ///< Suppresses growth or promotes quiescence
    INVASION,  ///< Modulates motility, ECM remodeling, EMT
    SECRETION, ///< Controls paracrine substrate secretion rates
    STRESS     ///< Oxygen/drug/oxidative stress response
};

// ----------------------------------------------------------------------------
// 3. GENE METADATA STRUCT
//    Encapsulates biological and simulation-level metadata for each gene node.
//    Stored in the compile-time GENE_INFO array indexed by GeneIndex.
// ----------------------------------------------------------------------------
struct GeneInfo
{
    const char*    name;                 ///< Human-readable gene symbol (use std::string_view for safe access)
    CellType       cell_type;            ///< Which cell type(s) express this gene
    FunctionalAxis functional_axis;      ///< Phenotype axis this gene modulates
    bool           is_druggable;         ///< Whether a clinical/experimental drug targets this node
    double         default_tumor_state;  ///< Initial Boolean/continuous state in tumor_cell (0.0–1.0)
    double         default_stroma_state; ///< Initial Boolean/continuous state in stromal_cell (0.0–1.0)
};

// ----------------------------------------------------------------------------
// 4. STATIC GENE METADATA ARRAY
//    Indexed by GeneIndex enum value.
//    Populated using C++17 aggregate initialization.
//
//    default_stroma_state is set to 0.0 for TUMOR-only genes (unused);
//    default_tumor_state  is set to 0.0 for STROMA-only genes (unused).
// ----------------------------------------------------------------------------
inline constexpr std::array<GeneInfo, GENE_COUNT> GENE_INFO = {{

    // Index 0 — KRAS
    // Constitutively active in virtually all PDAC (G12D/V hotspot mutation).
    // Drives cell cycle entry, metabolic reprogramming, anti-apoptosis.
    // Not directly druggable in this model (KRAS G12C inhibitors not canonical for PDAC).
    { "KRAS",   CellType::TUMOR,  FunctionalAxis::GROWTH,    false, 1.0, 0.0 },

    // Index 1 — MYC
    // Amplified/overexpressed downstream of KRAS. Intermediate default (0.5)
    // reflects heterogeneous expression; updated by Boolean network.
    { "MYC",    CellType::TUMOR,  FunctionalAxis::GROWTH,    false, 0.5, 0.0 },

    // Index 2 — CCND1
    // Cyclin D1; pairs with CDK4/6 to drive G1/S transition.
    // Freed from p16/CDKN2A inhibition in PDAC (CDKN2A LOF ~90%).
    // Druggable: palbociclib/ribociclib/abemaciclib (CDK4/6 inhibitors).
    // Initial state 0.5 (dynamically computed; rises quickly under KRAS/MYC drive).
    { "CCND1",  CellType::TUMOR,  FunctionalAxis::GROWTH,    true,  0.5, 0.0 },

    // Index 3 — TP53
    // Tumor suppressor; guardian of genome. Lost (mutated) in ~75% of PDAC.
    // When 0: apoptosis threshold elevated, DNA damage response blunted.
    { "TP53",   CellType::TUMOR,  FunctionalAxis::DEATH,     false, 0.0, 0.0 },

    // Index 4 — BCL_XL
    // Anti-apoptotic BCL-2 family member. Elevated in PDAC; navitoclax target.
    // High default (0.8) reflects KRAS-driven survival signaling.
    { "BCL_XL", CellType::TUMOR,  FunctionalAxis::DEATH,     true,  0.8, 0.0 },

    // Index 5 — SNAI1
    // Snail EMT transcription factor; fast (hours) CDH1 repressor.
    // First responder to TGF-beta signal; acts upstream of ZEB1.
    // Not directly druggable (no approved SNAI1 inhibitors).
    // Initial state 0.0 (quiescent epithelial; signal-driven).
    { "SNAI1",  CellType::TUMOR,  FunctionalAxis::INVASION,  false, 0.0, 0.0 },

    // Index 6 — CDKN2A
    // p16/INK4A; CDK4/6 inhibitor. Deleted in ~95% PDAC via homozygous deletion.
    // When 0: G1/S checkpoint fails, proliferation rate uninhibited.
    { "CDKN2A", CellType::TUMOR,  FunctionalAxis::BRAKING,   false, 0.0, 0.0 },

    // Index 7 — SMAD4
    // Canonical TGF-beta nuclear effector. Lost in ~55% PDAC.
    // When 0: TGF-beta cannot suppress proliferation; stroma crosstalk uncoupled.
    { "SMAD4",  CellType::TUMOR,  FunctionalAxis::BRAKING,   false, 0.0, 0.0 },

    // Index 8 — RB1
    // Retinoblastoma protein; E2F repressor at G1/S. Intermediate default (0.5)
    // reflects functional but partially bypassed state (CDKN2A loss).
    { "RB1",    CellType::TUMOR,  FunctionalAxis::BRAKING,   false, 0.5, 0.0 },

    // Index 9 — ZEB1
    // EMT transcription factor. Activated by TGF-beta/KRAS; represses CDH1.
    // Initially 0 (epithelial state); rises with TGF-beta exposure.
    { "ZEB1",   CellType::TUMOR,  FunctionalAxis::INVASION,  false, 0.0, 0.0 },

    // Index 10 — CDH1
    // E-cadherin; epithelial adhesion molecule. Repressed by ZEB1 during EMT.
    // Initially 1 (epithelial); decreases as ZEB1 rises.
    { "CDH1",   CellType::TUMOR,  FunctionalAxis::INVASION,  false, 1.0, 0.0 },

    // Index 11 — MMP2
    // Gelatinase B; secreted by both tumor and activated stroma.
    // Degrades type IV collagen and ECM; promotes invasion.
    // Druggable (MMP inhibitors in research; not yet clinical standard).
    { "MMP2",   CellType::BOTH,   FunctionalAxis::INVASION,  true,  0.2, 0.1 },

    // Index 12 — ACTA2
    // Alpha-smooth muscle actin; definitive CAF activation marker.
    // Quiescent PSCs express 0; activated CAFs express 1.
    // Stroma-only; no expression in tumor cells.
    { "ACTA2",  CellType::STROMA, FunctionalAxis::INVASION,  false, 0.0, 0.0 },

    // Index 13 — TGFB1
    // TGF-beta1 gene expression proxy (not the substrate field).
    // Expressed by both tumor (paracrine) and activated stroma (autocrine).
    // Druggable (galunisertib TGFbR1 inhibitor; fresolimumab anti-TGFb).
    { "TGFB1",  CellType::BOTH,   FunctionalAxis::SECRETION, true,  0.3, 0.1 },

    // Index 14 — SHH
    // Sonic Hedgehog ligand; secreted by tumor to activate stroma via SMO/GLI1.
    // Druggable: vismodegib/sonidegib (SMO inhibitors) block downstream signaling.
    { "SHH",    CellType::TUMOR,  FunctionalAxis::SECRETION, true,  0.3, 0.0 },

    // Index 15 — GLI1
    // Hedgehog pathway transcription factor; activated in stroma by SHH.
    // Drives ACTA2, HAS2, COL1A1 expression. Druggable (GANT61).
    // Stroma-only; initialized 0 (quiescent).
    { "GLI1",   CellType::STROMA, FunctionalAxis::SECRETION, true,  0.0, 0.0 },

    // Index 16 — HAS2
    // Hyaluronan synthase 2; synthesizes hyaluronic acid ECM component.
    // Driven by GLI1 in activated stroma; creates physical barrier to drug delivery.
    // Druggable: PEGPH20 (hyaluronidase) depletes HA to improve drug penetration.
    { "HAS2",   CellType::STROMA, FunctionalAxis::SECRETION, true,  0.0, 0.0 },

    // Index 17 — COL1A1
    // Collagen I alpha-1; major structural ECM protein deposited by CAFs.
    // Creates desmoplastic barrier; not directly druggable in this model.
    // Stroma-only; initialized 0 (quiescent).
    { "COL1A1", CellType::STROMA, FunctionalAxis::SECRETION, false, 0.0, 0.0 },

    // Index 18 — HIF1A
    // HIF-1alpha; master regulator of hypoxic response.
    // Activated in both cell types when oxygen < HYPOXIA_O2_THRESHOLD.
    // Drives anaerobic metabolism, VEGF, and drug resistance.
    { "HIF1A",  CellType::BOTH,   FunctionalAxis::STRESS,    false, 0.0, 0.0 },

    // Index 19 — NRF2
    // NRF2 (NFE2L2); antioxidant response transcription factor.
    // Activated by oxidative stress; confers drug resistance when high.
    // Tumor-only in this model; stroma NRF2 not explicitly modeled.
    { "NRF2",   CellType::TUMOR,  FunctionalAxis::STRESS,    false, 0.0, 0.0 },

    // Index 20 — ABCB1
    // MDR1/P-glycoprotein ATP-dependent drug efflux pump.
    // Upregulated by HIF1A and NRF2; reduces intracellular drug concentration.
    // Druggable (elacridar, tariquidar MDR1 inhibitors).
    { "ABCB1",  CellType::TUMOR,  FunctionalAxis::STRESS,    true,  0.0, 0.0 },

}};

// ----------------------------------------------------------------------------
// 5. COMPILE-TIME FALLBACK CONSTANTS
//
//    WARNING: These constants are FALLBACK DEFAULTS ONLY.
//    The simulation reads the authoritative values from PhysiCell_settings.xml
//    user_parameters at startup (via ThresholdConfig::load_from_xml below).
//    These fallbacks are used only if load_from_xml() is not called (e.g.,
//    unit tests, standalone tools).
//
//    KEEP THESE IN SYNC WITH PhysiCell_settings.xml user_parameters:
//      oxygen_hypoxia_threshold  = 15.0  (HIF1A activation, moderate hypoxia)
//      oxygen_necrosis_threshold =  5.0  (necrosis, severe hypoxia)
//      tgfb_activation_threshold =  0.1
//      shh_activation_threshold  =  0.05
// ----------------------------------------------------------------------------

/// O2 (mmHg) below which HIF1A activates — moderate hypoxia onset.
/// PHD2 Km for O2 is ~5 mmHg, but HIF1A accumulation begins ~15 mmHg in vivo.
/// Matches XML: oxygen_hypoxia_threshold = 15.0
inline constexpr double HYPOXIA_HIF1A_THRESHOLD = 15.0;

/// O2 (mmHg) below which necrosis rate is upregulated — severe / anoxic.
/// Matches XML: oxygen_necrosis_threshold = 5.0
inline constexpr double HYPOXIA_NECROSIS_THRESHOLD = 5.0;

/// tgfb concentration above which stromal CAF activation triggers.
/// Matches XML: tgfb_activation_threshold = 0.1
inline constexpr double TGFB_ACTIVATION_THRESHOLD = 0.1;

/// shh concentration above which GLI1 activates in stroma.
/// Matches XML: shh_activation_threshold = 0.05
inline constexpr double SHH_ACTIVATION_THRESHOLD = 0.05;

/// ecm_density above which stroma is considered dense/desmoplastic.
/// No XML counterpart; used only in tumor_cell.cpp motility reduction.
inline constexpr double ECM_HIGH_THRESHOLD = 0.7;

/// Neighbor contact fraction above which contact inhibition suppresses division.
inline constexpr double CONTACT_INHIBITION_THRESHOLD = 0.8;

/// ZEB1 level above which the cell is classified mesenchymal.
inline constexpr double EMT_ZEB1_THRESHOLD = 0.5;

/// NRF2 level above which drug resistance (ABCB1 upregulation) kicks in.
inline constexpr double STRESS_NRF2_THRESHOLD = 0.5;

// ----------------------------------------------------------------------------
// 6. ThresholdConfig — runtime threshold values loaded from XML
//
//    Single source of truth for all microenvironment thresholds.
//    Populated once at simulation setup via load_from_xml(), then passed
//    by const-ref to BooleanNetwork::update() and phenotype functions.
//
//    Fallback values (used when load_from_xml() is not called) match both
//    the compile-time constants above and PhysiCell_settings.xml.
//
//    Usage in custom.cpp setup:
//      ThresholdConfig cfg = ThresholdConfig::load_from_xml();
//      // store cfg globally or in a simulation context object
//
//    Usage in phenotype update:
//      bn.update(dt, o2, tgfb, shh, drug, cfg);
// ----------------------------------------------------------------------------
struct ThresholdConfig
{
    double hypoxia_hif1a_threshold   = HYPOXIA_HIF1A_THRESHOLD;   ///< O2 for HIF1A (mmHg)
    double hypoxia_necrosis_threshold= HYPOXIA_NECROSIS_THRESHOLD; ///< O2 for necrosis (mmHg)
    double tgfb_activation_threshold = TGFB_ACTIVATION_THRESHOLD;  ///< tgfb for CAF activation
    double shh_activation_threshold  = SHH_ACTIVATION_THRESHOLD;   ///< shh for GLI1 activation
    double ecm_high_threshold        = ECM_HIGH_THRESHOLD;          ///< dense ECM cutoff
    double contact_inhibition_threshold = CONTACT_INHIBITION_THRESHOLD;
    double emt_zeb1_threshold        = EMT_ZEB1_THRESHOLD;
    double stress_nrf2_threshold     = STRESS_NRF2_THRESHOLD;

    /// Construct with default values (compile-time fallbacks).
    ThresholdConfig() = default;

    /// Load authoritative values from PhysiCell XML user_parameters.
    /// Call once after PhysiCell has parsed the XML (after setup_microenvironment
    /// or at the start of setup_tissue).
    /// Returns a ThresholdConfig populated from parameters.doubles(); any
    /// parameter not present in the XML retains its default value.
    static ThresholdConfig load_from_xml();
};

#endif // GENE_DEFINITIONS_H
