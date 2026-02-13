#ifndef TEA_STATE_H
#define TEA_STATE_H

#include <ttak/math/bigreal.h>
#include <ttak/shared/shared.h>
#include <ttak/mem/owner.h>

typedef enum {
    VESSEL_BRITISH_TEAPOT,
    VESSEL_ZISHA,
    VESSEL_WHITE_PORCELAIN,
    VESSEL_CELADON,
    VESSEL_JAPANESE_CERAMICS,
    VESSEL_GLASS,
    VESSEL_GAIWAN_WHITE_PORCELAIN,
    VESSEL_GAIWAN_CELADON
} vessel_type_t;

typedef enum {
    BOIL_ELECTRIC,
    BOIL_POT_IRON,
    BOIL_POT_BRONZE,
    BOIL_POT_CLAY
} boil_method_t;

typedef enum {
    TEA_GREEN_NORMAL,
    TEA_WHITE,
    TEA_BLACK,
    TEA_OOLONG,
    TEA_YELLOW,
    TEA_PUERH,
    TEA_GREEN_GYOKURO,
    TEA_GREEN_SENCHA,
    TEA_GREEN_FUKAMUSHI,
    TEA_TIBETAN
} tea_type_t;

typedef enum {
    UNIT_METRIC,
    UNIT_US,
    UNIT_UK
} unit_system_t;

typedef enum {
    HARDNESS_TDS_PPM,
    HARDNESS_FH,
    HARDNESS_DH
} hardness_unit_t;

typedef struct {
    uint64_t start_time;
    ttak_bigreal_t current_temp;     // Internal always Celsius
    ttak_bigreal_t extraction_level; // Arbitrary units (e.g., 0.0 to 1.0)
    
    // Detailed Components (mg/g or relative) - Current Infusion
    ttak_bigreal_t ext_catechin;
    ttak_bigreal_t ext_amino_acid;
    ttak_bigreal_t ext_caffeine;
    ttak_bigreal_t ext_pectin;
    ttak_bigreal_t ext_polysaccharide;
    ttak_bigreal_t ext_lignin;       // Structural compounds released on rupture
    ttak_bigreal_t ext_aroma;        // Volatile aromatics captured during infusion

    // Accumulated from previous infusions
    ttak_bigreal_t accum_ext_catechin;
    ttak_bigreal_t accum_ext_amino_acid;
    ttak_bigreal_t accum_ext_caffeine;
    ttak_bigreal_t accum_ext_pectin;
    ttak_bigreal_t accum_ext_polysaccharide;
    ttak_bigreal_t accum_ext_lignin;
    ttak_bigreal_t accum_ext_aroma;
    
    // New Physics Model State
    ttak_bigreal_t structural_integrity; // 1.0 (intact) to 0.0 (destroyed)
    ttak_bigreal_t hydration_state;  // 0.0 to 1.0
    ttak_bigreal_t unfurling_state;  // 0.0 to 1.0 (leaf opening)
    ttak_bigreal_t leaf_density;     // g/cm^3
    ttak_bigreal_t astringency_rate; // d(astringency)/dt
    ttak_bigreal_t sweetness_rate;   // d(sweetness)/dt
    ttak_bigreal_t saturation_index; // Total saturation
    _Bool stop_signal;               // Stop trigger
    ttak_bigreal_t aroma_extraction_axis;  // mg/min of captured aromatics
    ttak_bigreal_t aroma_volatility_axis;  // mg/min being lost to volatility
    ttak_bigreal_t amino_depth_axis;       // Normalized amino body
    ttak_bigreal_t amino_vibrancy_axis;    // mg/min change in amino acids
    ttak_bigreal_t clarity_index;          // mg/ml of heavy colloids for style clamp
    
    // Unit Preferences
    unit_system_t unit_system;
    _Bool use_fahrenheit; // Fahrenheit for temp
    hardness_unit_t hardness_unit;
    
    // Derivatives for display
    ttak_bigreal_t ext_velocity;     // dC/dt
    ttak_bigreal_t bitter_accel;     // d2B/dt2

    // Tea Distributor Algorithm state
    ttak_bigreal_t distributor_potential;    // Remaining mg/ml that can still be extracted
    ttak_bigreal_t distributor_drag_ratio;   // Actual vs expected flow for rebalancing
    ttak_bigreal_t distributor_target_hint;  // Suggested allocation for the active infusion
    
    // Configuration
    tea_type_t tea_type;
    vessel_type_t vessel;
    boil_method_t boil_method;
    ttak_bigreal_t water_hardness; // PPM
    ttak_bigreal_t leaf_width;     // mm
    ttak_bigreal_t leaf_height;    // mm
    ttak_bigreal_t leaf_mass;      // grams
    ttak_bigreal_t water_volume_ml; // ml
    ttak_bigreal_t altitude_m;      // Altitude in meters
    uint32_t vintage_years;        // Age of tea
    _Bool is_ripe;                 // Pu-erh: Shu (Ripe) vs Sheng (Raw)
    _Bool has_butter;              // Tibetan tea butter layer
    _Bool boiling_in_pot;          // Whether tea is being boiled in a pot
    uint8_t heat_level;            // Pot heat (1-10)
    
    // Water Quality
    _Bool is_bottled;              // Mineral/Bottled water
    ttak_bigreal_t mineral_content[5]; // Ca, Mg, Na, K, etc.
    ttak_bigreal_t tds;            // Total Dissolved Solids
    _Bool is_limestone_bedrock;    // For tap water without TDS info
    
    // Infusion Tracking
    uint32_t num_infusions;
    uint32_t current_infusion;
    uint64_t cycle_start_time;
    ttak_bigreal_t target_extraction_for_cycle;
    ttak_bigreal_t cycle_start_extraction;
    _Bool cycle_active;

    // Status
    uint64_t last_update_time;
    _Bool temp_preserved;
    _Bool is_exhausted;
} tea_state_t;

TTAK_SHARED_DEFINE_WRAPPER(tea_state, tea_state_t)

#endif
