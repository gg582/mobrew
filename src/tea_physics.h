#ifndef TEA_PHYSICS_H
#define TEA_PHYSICS_H

#include <ttak/math/bigreal.h>
#include "tea_state.h"

// Constants for materials (simplified for simulation)
// These would ideally be looked up or configured.

typedef struct {
    ttak_bigreal_t specific_heat;        // J/g*K
    ttak_bigreal_t porosity;             // 0.0 - 1.0 factor affecting cooling/extraction
    ttak_bigreal_t thermal_conductivity; // W/m*K
} vessel_props_t;

typedef struct {
    double base_velocity;    // Base K for extraction
    double astringency_acc;  // Base acceleration factor
    double leaf_density;     // Typical density for the type
} tea_profile_t;

typedef struct {
    double s_catechin;      // Max soluble mg/g
    double s_amino_acid;
    double s_caffeine;
    double s_pectin;
    double s_polysaccharide;
    double s_aroma;

    // Extraction rate constants at 100C (min^-1)
    double k_catechin;
    double k_amino_acid;
    double k_caffeine;
    double k_pectin;
    double k_polysaccharide;
    double k_aroma;

    // Temperature sensitivity (Activation Energy factor / R)
    // Higher means more sensitive to temperature
    double sens_catechin;
    double sens_amino_acid;
    double sens_caffeine;
    double sens_pectin;
    double sens_polysaccharide;
    double sens_aroma;

    // Base volatility loss multiplier for aromatics
    double aroma_volatility_base;
} tea_comp_profile_t;

void get_vessel_props(vessel_type_t type, vessel_props_t *props);
void get_tea_profile(tea_type_t type, tea_profile_t *profile);
void get_tea_comp_profile(tea_type_t type, tea_comp_profile_t *profile);

// Main physics step function
// Advances state from last_update_time to now.
void physics_simulate_step(tea_state_t *state, uint64_t now);

// New helper to simulate lag in discrete steps
void physics_simulate_lag(tea_state_t *state, uint64_t lag_ms);

// Backtrack functions for initial setup
void physics_backtrack_cooling(ttak_bigreal_t *temp_out, ttak_bigreal_t *start_temp, uint64_t duration_ms);
void physics_backtrack_extraction(ttak_bigreal_t *extract_out, const ttak_bigreal_t *temp_profile, uint64_t duration_ms);

// Helper to get extraction rate at a specific state (for derivative display)
void physics_get_extraction_rate(ttak_bigreal_t *rate_out, const tea_state_t *state);

// Helper to get bitterness acceleration
void physics_get_bitterness_accel(ttak_bigreal_t *accel_out, const tea_state_t *state);

// Tea Distributor algorithm helper
double tea_distributor_adjust_target(tea_state_t *state, double base_slice_pct);

// Cleanup helper
void tea_state_init(tea_state_t *state);
void tea_state_cleanup(void *ptr);

// Conversion helpers
void set_br_d(ttak_bigreal_t *out, double val);
double get_double_br(const ttak_bigreal_t *br);

#endif
