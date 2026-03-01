#ifndef CUSTOM_H
#define CUSTOM_H

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"

#include "stromal_cell.h"
#include "tumor_cell.h"

#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Global substrate indices (set in setup_microenvironment()).
// Exposed so phenotype modules can reuse a single source of truth.
// ---------------------------------------------------------------------------
extern int oxygen_index;
extern int tgfb_index;
extern int shh_index;
extern int ecm_index;
extern int drug_index;

struct ECMInterventionState
{
    bool ha_degrade_active = false;
    bool col_degrade_active = false;
    double ha_degrade_strength = 0.0;
    double col_degrade_strength = 0.0;
};

extern ECMInterventionState intervention_state;

// ---------------------------------------------------------------------------
// Phase-gated module execution
// ---------------------------------------------------------------------------
enum class ModulePhase
{
    SENSING,
    DECISION,
    WRITE
};

struct ModuleTraceEntry
{
    std::string module_name;
    ModulePhase phase;
    std::vector<std::string> field_reads;
    std::vector<std::string> field_writes;
};

// ---------------------------------------------------------------------------
// PhysiCell custom module entry points
// ---------------------------------------------------------------------------
void setup_microenvironment(void);
void create_cell_types(void);
void setup_tissue(void);

// Dispatcher for phenotype update.
void custom_function(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt);

// SENSING PHASE
void module1_oxygen_sensing(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt, ModulePhase phase);
void module3_stromal_activation(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt, ModulePhase phase);
void module7_drug_response(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt, ModulePhase phase);

// DECISION PHASE
void module5_emt_engine(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt, ModulePhase phase);
void module4_proliferation_death(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt, ModulePhase phase);

// WRITE PHASE
void module2_paracrine_secretion(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt, ModulePhase phase);
void module6_ecm_production(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt, ModulePhase phase);
void module8_mechanical_compaction(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt, ModulePhase phase);

// Optional coloring callback used by standard PhysiCell main.cpp templates.
std::vector<std::string> my_coloring_function(PhysiCell::Cell* pCell);

// ---------------------------------------------------------------------------
// ECM barrier microenvironment modifier
// ---------------------------------------------------------------------------
void ecm_dependent_diffusion(double dt);
void ecm_dependent_diffusion_solver(BioFVM::Microenvironment& M, double dt);
void register_ecm_dependent_diffusion_solver(void);
void update_ecm_effective_diffusion_coefficients(BioFVM::Microenvironment& M);
double get_effective_diffusion_coefficient(int substrate_index, int voxel_index);
double get_ecm_collagen_fraction(int voxel_index);
double get_local_mechanical_stiffness(int voxel_index);

// Voxel-level ECM composition storage used by Module 6.
void reset_ecm_ha_fraction_field(double default_fraction);
double get_ecm_ha_fraction(int voxel_index);
void set_ecm_ha_fraction(int voxel_index, double value);

// Module execution tracing used by order/phase verification tests.
void set_module_trace_enabled(bool enabled);
void clear_module_trace_log(void);
std::vector<ModuleTraceEntry> get_module_trace_log(void);

#endif // CUSTOM_H
