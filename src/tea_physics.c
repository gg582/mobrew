#include "tea_physics.h"
#include <ttak/math/calculus.h>
#include <ttak/math/bigreal.h>
#include <ttak/timing/timing.h>
#include <math.h> 
#include <stdio.h>

// Helper to init bigreal from double (approx)
void set_br_d(ttak_bigreal_t *out, double val) {
    int64_t i = (int64_t)(val * 1000000.0);
    ttak_bigint_set_u64(&out->mantissa, (i < 0 ? -i : i), ttak_get_tick_count());
    out->exponent = -6; 
    out->mantissa.is_negative = (val < 0);
}

// Helper to export double from bigreal (approx)
double get_double_br(const ttak_bigreal_t *br) {
    uint64_t val = 0;
    // This is a hacky way since we know internal structure from previous steps
    // or use export if available. But export_u64 gives integer part.
    // For this prototype, we rely on the fact that we constructed them with exponent -6 often.
    // If exponent is 0, it is integer.
    
    // Check if we can use export_u64
    if (ttak_bigint_export_u64(&br->mantissa, &val)) {
        double d = (double)val;
        if (br->mantissa.is_negative) d = -d; // Handle sign
        if (br->exponent == -6) return d / 1000000.0;
        if (br->exponent == 0) return d;
        return d * pow(10, br->exponent);
    }
    return 0.0;
}

void get_vessel_props(vessel_type_t type, vessel_props_t *props) {
    if (!props) return;
    *props = (vessel_props_t){0};
    double sh = 0.8; 
    double por = 0.1; 
    double tc = 1.0; // W/m*K default
    
    switch (type) {
        case VESSEL_BRITISH_TEAPOT: sh = 0.84; por = 0.05; tc = 1.1; break;
        case VESSEL_ZISHA:          sh = 0.90; por = 0.20; tc = 2.5; break;
        case VESSEL_WHITE_PORCELAIN: sh = 0.80; por = 0.01; tc = 1.5; break;
        case VESSEL_CELADON:        sh = 0.82; por = 0.02; tc = 1.6; break;
        case VESSEL_JAPANESE_CERAMICS: sh = 0.85; por = 0.10; tc = 1.8; break;
        case VESSEL_GLASS:           sh = 0.75; por = 0.00; tc = 0.9; break;
        case VESSEL_GAIWAN_WHITE_PORCELAIN: sh = 0.78; por = 0.01; tc = 1.7; break;
        case VESSEL_GAIWAN_CELADON:        sh = 0.80; por = 0.02; tc = 1.8; break;
    }
    
    set_br_d(&props->specific_heat, sh);
    set_br_d(&props->porosity, por);
    set_br_d(&props->thermal_conductivity, tc);
}

void get_tea_profile(tea_type_t type, tea_profile_t *profile) {
    // Default density: 0.5 g/cm^3
    profile->leaf_density = 0.5;
    
    switch (type) {
        case TEA_GREEN_NORMAL:    profile->base_velocity = 0.120; profile->astringency_acc = 0.002; profile->leaf_density = 0.45; break;
        case TEA_WHITE:           profile->base_velocity = 0.080; profile->astringency_acc = 0.001; profile->leaf_density = 0.35; break;
        case TEA_BLACK:           profile->base_velocity = 0.160; profile->astringency_acc = 0.005; profile->leaf_density = 0.55; break;
        case TEA_OOLONG:          profile->base_velocity = 0.110; profile->astringency_acc = 0.003; profile->leaf_density = 0.60; break;
        case TEA_YELLOW:          profile->base_velocity = 0.100; profile->astringency_acc = 0.002; profile->leaf_density = 0.48; break;
        case TEA_PUERH:           profile->base_velocity = 0.180; profile->astringency_acc = 0.001; profile->leaf_density = 0.70; break;
        case TEA_GREEN_GYOKURO:   profile->base_velocity = 0.140; profile->astringency_acc = 0.004; profile->leaf_density = 0.40; break;
        case TEA_GREEN_SENCHA:    profile->base_velocity = 0.130; profile->astringency_acc = 0.003; profile->leaf_density = 0.42; break;
        case TEA_GREEN_FUKAMUSHI: profile->base_velocity = 0.170; profile->astringency_acc = 0.004; profile->leaf_density = 0.38; break;
        case TEA_TIBETAN:         profile->base_velocity = 0.060; profile->astringency_acc = 0.000; profile->leaf_density = 0.80; break;
        default:                  profile->base_velocity = 0.100; profile->astringency_acc = 0.002; break;
    }
}

void get_tea_comp_profile(tea_type_t type, tea_comp_profile_t *profile) {
    // Increased constants (5-10x) to reflect 1-3 minute infusions
    profile->s_catechin = 120.0;
    profile->s_amino_acid = 20.0;
    profile->s_caffeine = 35.0;
    profile->s_pectin = 18.0;
    profile->s_polysaccharide = 22.0;
    profile->s_aroma = 18.0;

    profile->k_catechin = 0.40;
    profile->k_amino_acid = 0.70;
    profile->k_caffeine = 0.50;
    profile->k_pectin = 0.30;
    profile->k_polysaccharide = 0.25;
    profile->k_aroma = 0.55;

    profile->sens_catechin = 3.0;      
    profile->sens_amino_acid = 0.8;   
    profile->sens_caffeine = 2.0;
    profile->sens_pectin = 2.5;
    profile->sens_polysaccharide = 1.5;
    profile->sens_aroma = 1.2;
    profile->aroma_volatility_base = 0.30;

    switch (type) {
        case TEA_GREEN_NORMAL:
        case TEA_GREEN_SENCHA:
        case TEA_GREEN_FUKAMUSHI:
            profile->s_catechin = 150.0;
            profile->s_amino_acid = 30.0;
            profile->s_caffeine = 30.0;
            profile->s_aroma = 22.0;
            // High amino acid extraction even at low temp
            profile->k_amino_acid = 0.90; 
            profile->sens_amino_acid = 0.5; // Very low sensitivity, extracts well at cold
            // Catechins still temp sensitive but less than black
            profile->sens_catechin = 2.5;
            profile->k_aroma = 0.75;
            profile->aroma_volatility_base = 0.22;
            break;
        case TEA_GREEN_GYOKURO:
            profile->s_catechin = 130.0;
            profile->s_amino_acid = 60.0; 
            profile->k_amino_acid = 1.20; // Very fast amino acid release
            profile->sens_amino_acid = 0.3; // Extracts heavily at 50-60C
            profile->sens_catechin = 3.5;   // Suppress bitterness at low temp
            profile->s_aroma = 28.0;
            profile->k_aroma = 0.90;
            profile->aroma_volatility_base = 0.18;
            break;
        case TEA_WHITE:
            profile->s_catechin = 140.0;
            profile->s_amino_acid = 35.0;
            profile->s_aroma = 24.0;
            // Aggressive low temp extraction
            profile->k_catechin = 0.50; 
            profile->sens_catechin = 1.8; // Reduced temp sensitivity
            profile->sens_amino_acid = 0.6;
            profile->k_aroma = 0.65;
            profile->aroma_volatility_base = 0.20;
            break;
        case TEA_YELLOW:
            profile->s_catechin = 135.0;
            profile->s_amino_acid = 28.0;
            profile->s_aroma = 20.0;
            // Similar to green but slightly more oxidized/menstruated
            profile->k_catechin = 0.45;
            profile->sens_catechin = 2.0;
            profile->sens_amino_acid = 0.7;
            profile->k_aroma = 0.60;
            break;
        case TEA_BLACK:
            profile->s_catechin = 100.0; 
            profile->s_amino_acid = 10.0;
            profile->s_caffeine = 45.0;
            profile->s_aroma = 16.0;
            profile->k_catechin = 0.60;
            profile->k_caffeine = 0.60;
            profile->k_aroma = 0.50;
            profile->aroma_volatility_base = 0.38;
            break;
        case TEA_PUERH:
            profile->s_catechin = 90.0;
            profile->s_polysaccharide = 40.0; 
            profile->k_polysaccharide = 0.45;
            profile->s_aroma = 12.0;
            profile->k_aroma = 0.35;
            profile->aroma_volatility_base = 0.42;
            break;
        default: break;
    }
}

static double get_aroma_perception_floor(tea_type_t type) {
    switch (type) {
        case TEA_PUERH:
        case TEA_TIBETAN:
            return 0.30;
        case TEA_BLACK:
        case TEA_OOLONG:
            return 0.40;
        case TEA_GREEN_GYOKURO:
        case TEA_WHITE:
            return 0.36;
        default:
            return 0.35;
    }
}

static double get_aroma_balance_tolerance(tea_type_t type) {
    switch (type) {
        case TEA_PUERH:
        case TEA_TIBETAN:
            return 0.08;
        case TEA_BLACK:
        case TEA_OOLONG:
            return 0.05;
        default:
            return 0.04;
    }
}

static double get_aroma_signal_floor(tea_type_t type) {
    switch (type) {
        case TEA_PUERH:
        case TEA_TIBETAN:
            return 0.45;
        case TEA_BLACK:
        case TEA_OOLONG:
            return 0.70;
        default:
            return 0.60;
    }
}

static double get_aroma_balance_min_time_ms(tea_type_t type) {
    switch (type) {
        case TEA_PUERH:
            return 25000.0;
        case TEA_TIBETAN:
            return 30000.0;
        case TEA_BLACK:
        case TEA_OOLONG:
            return 14000.0;
        default:
            return 10000.0;
    }
}

static double get_aroma_stop_flux_min_extraction(tea_type_t type) {
    (void)type;
    return 0.06;
}

static double get_aroma_stop_flux_min_volatility(tea_type_t type) {
    (void)type;
    return 0.04;
}

static double clampd(double val, double min_v, double max_v) {
    if (val < min_v) return min_v;
    if (val > max_v) return max_v;
    return val;
}

static double compute_ratio_diffusion_guard_ms(const tea_state_t *state) {
    if (!state) return 5000.0;

    double leaf_mass = get_double_br(&state->leaf_mass);
    if (leaf_mass <= 0.05) leaf_mass = 0.05;

    double water_volume = get_double_br(&state->water_volume_ml);
    if (water_volume <= 0.0) water_volume = 100.0;

    double ratio = water_volume / leaf_mass; // ml per gram

    const double gongfu_ratio = 18.0;
    const double western_ratio = 38.0;
    const double min_time_sec = 5.0;
    const double max_time_sec = 28.0;

    double normalized = (ratio - gongfu_ratio) / (western_ratio - gongfu_ratio);
    normalized = clampd(normalized, 0.0, 1.0);
    double base_time_sec = min_time_sec + normalized * (max_time_sec - min_time_sec);

    if (ratio > western_ratio) {
        double excess = clampd((ratio - western_ratio) / 20.0, 0.0, 1.0);
        base_time_sec += excess * 6.0;
    }

    double hydration = clampd(get_double_br(&state->hydration_state), 0.0, 1.0);
    double hydration_relief = 1.0 - 0.35 * hydration;
    base_time_sec *= clampd(hydration_relief, 0.65, 1.0);

    if (base_time_sec < 4.5) base_time_sec = 4.5;
    return base_time_sec * 1000.0;
}

static double compute_leaf_geometry_factor(const tea_state_t *state, const tea_profile_t *tp, double *length_factor_out, double *leaf_count_out) {
    double width = get_double_br(&state->leaf_width);
    if (width <= 0.0) width = 6.0; // mm
    double height = get_double_br(&state->leaf_height);
    if (height <= 0.0) height = 15.0;

    double dominant = (height > width) ? height : width;
    double length_factor = tanh(dominant / 45.0) + 0.5; // 0.5 - 1.5 range
    if (length_factor_out) *length_factor_out = length_factor;

    double area_cm2 = (width * height) / 100.0; // mm^2 to cm^2
    if (area_cm2 < 0.4) area_cm2 = 0.4;

    double density = get_double_br(&state->leaf_density);
    if (density <= 0.0) density = tp->leaf_density;
    if (density < 0.25) density = 0.25;

    const double thickness_cm = 0.03; // 0.3 mm wet thickness assumption
    double per_leaf_mass = density * area_cm2 * thickness_cm;
    if (per_leaf_mass < 0.03) per_leaf_mass = 0.03;

    double leaf_mass = get_double_br(&state->leaf_mass);
    if (leaf_mass <= 0.0) leaf_mass = 5.0;

    double leaf_count = leaf_mass / per_leaf_mass;
    if (leaf_count < 1.0) leaf_count = 1.0;
    if (leaf_count_out) *leaf_count_out = leaf_count;

    double geometry_factor = 0.6 + 0.25 * (leaf_count / (leaf_count + 4.0)) + 0.15 * length_factor;
    return clampd(geometry_factor, 0.6, 1.6);
}

static void tea_distributor_refresh_snapshot(tea_state_t *state, const tea_profile_t *tp, const tea_comp_profile_t *cp) {
    if (!state || !tp || !cp) return;

    double length_factor = 1.0;
    double geometry_factor = compute_leaf_geometry_factor(state, tp, &length_factor, NULL);

    double total_solubles_per_g = cp->s_catechin + cp->s_amino_acid + cp->s_caffeine + cp->s_pectin + cp->s_polysaccharide;
    double leaf_mass = get_double_br(&state->leaf_mass);
    if (leaf_mass <= 0.0) leaf_mass = 5.0;

    double total_solids = total_solubles_per_g * leaf_mass;
    double extracted_pct = get_double_br(&state->extraction_level) / 100.0;
    if (extracted_pct < 0.0) extracted_pct = 0.0;
    if (extracted_pct > 1.0) extracted_pct = 1.0;

    double remaining_solids = total_solids * (1.0 - extracted_pct);
    if (remaining_solids < 0.0) remaining_solids = 0.0;

    double water_volume = get_double_br(&state->water_volume_ml);
    if (water_volume < 30.0) water_volume = 30.0;

    double potential_conc = (remaining_solids * geometry_factor) / water_volume;
    set_br_d(&state->distributor_potential, potential_conc);

    double hydration = get_double_br(&state->hydration_state);
    if (hydration < 0.05) hydration = 0.05;
    double expected_rate = tp->base_velocity * (0.8 + 0.2 * geometry_factor) * sqrt(hydration);
    if (expected_rate < 0.001) expected_rate = 0.001;

    double ext_velocity_raw = get_double_br(&state->ext_velocity);
    double actual_rate = (ext_velocity_raw > 0.0) ? (ext_velocity_raw / 100.0) : expected_rate;
    double drag_ratio = actual_rate / expected_rate;
    drag_ratio = clampd(drag_ratio, 0.4, 1.8);
    set_br_d(&state->distributor_drag_ratio, drag_ratio);
}

// Internal helper to calculate current K (rate constant) based on state for a component
static double calculate_component_k(const tea_state_t *state, double base_k, double sens, double hydration_factor) {
    // 1. Temperature Correction using Arrhenius-like model
    // k = base_k * exp(-Ea / (R * T)) -> simplified to k = base_k * exp(sens * (T - 100) / 30)
    double temp = get_double_br(&state->current_temp);
    double k = base_k * exp(sens * (temp - 100.0) / 30.0);

    // 2. Hydration-Diffusion Coupling (MFR Model)
    // k(t) = k_base(T) * sqrt(H(t))
    k *= hydration_factor;
    
    // 3. Water Quality Correction
    double tds = get_double_br(&state->tds);
    if (tds > 100.0) {
        double water_factor = 1.0 - ((tds - 100.0) / 2000.0); 
        if (water_factor < 0.7) water_factor = 0.7;
        k *= water_factor;
    }
    
    return k;
}

void physics_simulate_step(tea_state_t *state, uint64_t now) {
    if (now <= state->last_update_time) return;
    
    double dt_min = (now - state->last_update_time) / 1000.0 / 60.0;
    if (dt_min <= 0) return;

    uint64_t elapsed_ms = now - state->cycle_start_time;

    // 0. Hydration and Unfurling State Management
    double h = get_double_br(&state->hydration_state);
    double u = get_double_br(&state->unfurling_state);
    
    tea_profile_t tp;
    get_tea_profile(state->tea_type, &tp);
    double temp = get_double_br(&state->current_temp);
    double density = get_double_br(&state->leaf_density);
    if (density <= 0) density = tp.leaf_density;

    // Advanced Hydration Logic
    double k_hyd = 0.5;
    double temp_factor_hyd = 1.0;
    
    // Unfurling Logic
    double k_unfurl = 1.0; 

    switch (state->tea_type) {
        case TEA_GREEN_NORMAL:
        case TEA_GREEN_SENCHA:
        case TEA_GREEN_FUKAMUSHI:
        case TEA_YELLOW:
        case TEA_WHITE:
        case TEA_GREEN_GYOKURO:
                // Tender leaves: Hydrate fast and unfurl fast even at low temp
                k_hyd = 1.2;
                k_unfurl = 2.0; 
                // Reduced temperature dependency for hydration (can hydrate at 40-50C)
                // Base floor of 0.15 ensures activity even when cold
                temp_factor_hyd = exp((temp - 100.0) / 45.0) + 0.15;
                break;
            case TEA_OOLONG:
            case TEA_BLACK:
                // Rolled/Twisted leaves: Slower hydration and unfurling
                k_hyd = 0.6;
                k_unfurl = 0.8;
                temp_factor_hyd = exp((temp - 100.0) / 20.0);
                break;
            case TEA_PUERH:
            case TEA_TIBETAN:
                // Compressed: Very slow initially
                k_hyd = 0.3;
                k_unfurl = 0.4;
                temp_factor_hyd = exp((temp - 100.0) / 15.0);
                break;    }

    // Hydrate
    double dh_dt = k_hyd * (1.0 - h) * temp_factor_hyd / density;
    h += dh_dt * dt_min;
    if (h > 1.0) h = 1.0;
    set_br_d(&state->hydration_state, h);

    // Unfurl (follows hydration)
    // Mechanical expansion driven by turgor pressure (hydration) and heat
    double du_dt = k_unfurl * (h - u) * (temp / 100.0);
    // Acceleration as it opens (avalanche effect)
    if (u > 0.2) du_dt *= 1.2;
    
    u += du_dt * dt_min;
    if (u > 1.0) u = 1.0;
    set_br_d(&state->unfurling_state, u);

    // Effective Surface Area Factor
    // u=0 -> 20% exposed (surface only), u=1 -> 100% exposed
    double surface_factor = 0.2 + 0.8 * u;
    
    // Combined Factor for Extraction
    double hydration_factor = sqrt(h) * surface_factor;

    // 1. Thermal Dynamics & Boiling Detection
    vessel_props_t vprops;
    get_vessel_props(state->vessel, &vprops);
    double tc = get_double_br(&vprops.thermal_conductivity);
    double sh = get_double_br(&vprops.specific_heat);
    
    double alt = get_double_br(&state->altitude_m);
    double boiling_point = 100.0 - (alt / 285.0);
    if (boiling_point < 70.0) boiling_point = 70.0;

    _Bool is_active_boiling = (temp >= boiling_point - 0.5 && state->boiling_in_pot && state->heat_level > 0);

    // Turbulence Factor: Forced Convection during Active Boiling
    // Destroys boundary layer, increasing mass transfer coefficient significantly.
    double turbulence_factor = 1.0;
    if (is_active_boiling) {
        turbulence_factor = 3.0; 
        
        // Vaporization & Concentration Model
        // Evaporation rate scales with heat level.
        // Simplified: Level 10 = ~5 ml/min evaporation
        double evap_rate = 0.5 * state->heat_level; // ml/min
        double vol_lost = evap_rate * dt_min;
        double current_vol = get_double_br(&state->water_volume_ml);
        
        if (current_vol > vol_lost) {
            double new_vol = current_vol - vol_lost;
            set_br_d(&state->water_volume_ml, new_vol);
        }
    }

    // Cooling Logic (applied if not actively maintaining boil)
    if (!is_active_boiling) {
        double k_cool = 0.018 * (tc / sh); 
        if (state->has_butter) k_cool *= 0.15;
        double new_temp = temp * (1.0 - k_cool * dt_min);
        if (new_temp < 20.0) new_temp = 20.0;
        set_br_d(&state->current_temp, new_temp);
    } else {
        // Clamp to boiling point
        set_br_d(&state->current_temp, boiling_point);
    }

    // 2. Cell-Wall Rupture Model (Decoction Specific)
    double integrity = get_double_br(&state->structural_integrity);
    if (is_active_boiling) {
        // Integrity degrades over time under boiling conditions
        // Rate depends on leaf density (lighter leaves break faster)
        double k_rupture = 0.15; // Base rupture rate (per min)
        integrity -= k_rupture * dt_min;
        if (integrity < 0.0) integrity = 0.0;
        set_br_d(&state->structural_integrity, integrity);
    }

    // 3. Precision Extraction Tracking
    tea_comp_profile_t cp;
    get_tea_comp_profile(state->tea_type, &cp);

    double e_cat_curr = get_double_br(&state->ext_catechin);
    double e_ami_curr = get_double_br(&state->ext_amino_acid);
    double e_caf_curr = get_double_br(&state->ext_caffeine);
    double e_pec_curr = get_double_br(&state->ext_pectin);
    double e_pol_curr = get_double_br(&state->ext_polysaccharide);
    double e_lig_curr = get_double_br(&state->ext_lignin);
    double e_aroma_curr = get_double_br(&state->ext_aroma);

    double e_cat = e_cat_curr + get_double_br(&state->accum_ext_catechin);
    double e_ami = e_ami_curr + get_double_br(&state->accum_ext_amino_acid);
    double e_caf = e_caf_curr + get_double_br(&state->accum_ext_caffeine);
    double e_pec = e_pec_curr + get_double_br(&state->accum_ext_pectin);
    double e_pol = e_pol_curr + get_double_br(&state->accum_ext_polysaccharide);
    double e_lig = e_lig_curr + get_double_br(&state->accum_ext_lignin);
    double e_aroma = e_aroma_curr + get_double_br(&state->accum_ext_aroma);

    double lag_penalty = 1.0;
    if (elapsed_ms < 30000 && !is_active_boiling) {
         if (elapsed_ms < 10000) lag_penalty = 0.05;
         else lag_penalty = 0.05 + 0.95 * ((double)(elapsed_ms - 10000) / 20000.0);
    }

    // Vessel-specific extraction characteristics
    double vessel_extraction_factor = 1.0;
    if (state->vessel == VESSEL_GAIWAN_WHITE_PORCELAIN || state->vessel == VESSEL_GAIWAN_CELADON) {
        // Gaiwan: Higher surface area/agitation efficiency
        vessel_extraction_factor = 1.25;
    }

    // Apply Turbulence Factor and Vessel Factor to K values
    double k_cat = calculate_component_k(state, cp.k_catechin, cp.sens_catechin, hydration_factor) * lag_penalty * turbulence_factor * vessel_extraction_factor;
    double k_ami = calculate_component_k(state, cp.k_amino_acid, cp.sens_amino_acid, hydration_factor) * lag_penalty * turbulence_factor * vessel_extraction_factor;
    double k_caf = calculate_component_k(state, cp.k_caffeine, cp.sens_caffeine, hydration_factor) * lag_penalty * turbulence_factor * vessel_extraction_factor;
    double k_pol = calculate_component_k(state, cp.k_polysaccharide, cp.sens_polysaccharide, hydration_factor) * lag_penalty * turbulence_factor * vessel_extraction_factor;
    double k_aroma = calculate_component_k(state, cp.k_aroma, cp.sens_aroma, hydration_factor * surface_factor) * lag_penalty * turbulence_factor * vessel_extraction_factor;
    
    double thermal_pressure = sh / 0.8; 
    double k_pec = calculate_component_k(state, cp.k_pectin * thermal_pressure, cp.sens_pectin, hydration_factor) * lag_penalty * turbulence_factor * vessel_extraction_factor;

    double de_aroma = k_aroma * (cp.s_aroma - e_aroma) * dt_min;
    if (de_aroma < 0.0) de_aroma = 0.0;

    double openness = 1.0;
    switch (state->vessel) {
        case VESSEL_GLASS:
        case VESSEL_BRITISH_TEAPOT:
            openness = 1.25; break;
        case VESSEL_GAIWAN_WHITE_PORCELAIN:
        case VESSEL_GAIWAN_CELADON:
            openness = 1.35; break;
        case VESSEL_ZISHA:
            openness = 1.10; break;
        default: break;
    }
    if (state->has_butter) openness *= 0.6;

    double volatility_rate = cp.aroma_volatility_base * openness;
    double vapor_drive = exp((temp - 65.0) / 22.0);
    if (vapor_drive < 0.2) vapor_drive = 0.2;
    if (is_active_boiling) vapor_drive *= 1.4;
    volatility_rate *= vapor_drive;

    double available_aroma = e_aroma_curr + de_aroma;
    if (available_aroma < 0.0) available_aroma = 0.0;
    double aroma_vol_loss = available_aroma * volatility_rate * dt_min;
    if (aroma_vol_loss > available_aroma) aroma_vol_loss = available_aroma;
    double aroma_retained = available_aroma - aroma_vol_loss;
    if (aroma_retained < 0.0) aroma_retained = 0.0;
    const double aroma_perception_floor_ratio = get_aroma_perception_floor(state->tea_type);
    double aroma_perception_floor_mass = 0.0;
    if (cp.s_aroma > 0.0 && aroma_perception_floor_ratio > 0.0) {
        aroma_perception_floor_mass = cp.s_aroma * aroma_perception_floor_ratio;
    }
    const _Bool aroma_floor_met = (cp.s_aroma <= 0.0) || 
                                  (aroma_retained + 1e-6 >= aroma_perception_floor_mass);

    // Lignin Extraction (Triggered by Rupture)
    double k_lig = 0.0;
    double s_lignin = 50.0; // Max available lignin
    if (integrity < 0.3) {
        // Internal Solutes Release
        k_lig = 0.2 * (1.0 - integrity / 0.3) * turbulence_factor; // Accelerates as integrity fails
    }

    double de_cat = k_cat * (cp.s_catechin - e_cat) * dt_min;
    double de_ami = k_ami * (cp.s_amino_acid - e_ami) * dt_min;
    double de_caf = k_caf * (cp.s_caffeine - e_caf) * dt_min;
    double clarity_bias = 1.0;
    if (!is_active_boiling) {
        double sat_prev = clampd(1.0 - get_double_br(&state->saturation_index), 0.15, 0.95);
        clarity_bias = 0.55 + 0.45 * sat_prev;
    }

    double de_pec = k_pec * (cp.s_pectin - e_pec) * dt_min * clarity_bias;
    double de_pol = k_pol * (cp.s_polysaccharide - e_pol) * dt_min * clarity_bias;
    double de_lig = k_lig * (s_lignin - e_lig) * dt_min;

    set_br_d(&state->ext_catechin, e_cat_curr + de_cat);
    set_br_d(&state->ext_amino_acid, e_ami_curr + de_ami);
    set_br_d(&state->ext_caffeine, e_caf_curr + de_caf);
    set_br_d(&state->ext_pectin, e_pec_curr + de_pec);
    set_br_d(&state->ext_polysaccharide, e_pol_curr + de_pol);
    set_br_d(&state->ext_lignin, e_lig_curr + de_lig);
    set_br_d(&state->ext_aroma, aroma_retained);

    double current_vol = get_double_br(&state->water_volume_ml);
    if (current_vol < 10.0) current_vol = 10.0;
    double l_mass = get_double_br(&state->leaf_mass);
    if (l_mass <= 0.0) l_mass = 5.0;
    double heavy_solutes_conc = ( (e_pec + de_pec) * 1.5 + (e_pol + de_pol) * 1.0 + (e_lig + de_lig) * 2.0 ) * l_mass / current_vol;
    set_br_d(&state->clarity_index, heavy_solutes_conc);

    double aroma_ext_rate = (dt_min > 0.0) ? de_aroma / dt_min : 0.0;
    double aroma_vol_rate = (dt_min > 0.0) ? aroma_vol_loss / dt_min : 0.0;
    set_br_d(&state->aroma_extraction_axis, aroma_ext_rate);
    set_br_d(&state->aroma_volatility_axis, aroma_vol_rate);

    double aroma_capture_ratio = 0.0;
    if (cp.s_aroma > 0.0) {
        aroma_capture_ratio = clampd(aroma_retained / cp.s_aroma, 0.0, 1.2);
    }

    double aroma_guard_margin = aroma_ext_rate - aroma_vol_rate;
    double aroma_balance_ratio = 0.0;
    if (aroma_ext_rate > 1e-6) {
        aroma_balance_ratio = aroma_vol_rate / (aroma_ext_rate + 1e-6);
    }
    double min_flux_ext = get_aroma_stop_flux_min_extraction(state->tea_type);
    double min_flux_vol = get_aroma_stop_flux_min_volatility(state->tea_type);
    _Bool live_flux_gate = (aroma_ext_rate >= min_flux_ext && aroma_vol_rate >= min_flux_vol);

    double projected_aroma_ext_rate = aroma_ext_rate;
    double projected_aroma_vol_rate = aroma_vol_rate;
    double projected_aroma_balance = aroma_balance_ratio;
    double projected_guard_margin = aroma_ext_rate - aroma_vol_rate;
    if (dt_min > 0.0) {
        double remaining_aroma = cp.s_aroma - (e_aroma + de_aroma);
        if (remaining_aroma < 0.0) remaining_aroma = 0.0;
        double projected_de_aroma = k_aroma * remaining_aroma * dt_min;
        if (projected_de_aroma < 0.0) projected_de_aroma = 0.0;

        double next_ext_rate = projected_de_aroma / dt_min;
        if (next_ext_rate < 0.0) next_ext_rate = 0.0;

        double next_available_aroma = aroma_retained + projected_de_aroma;
        if (next_available_aroma < 0.0) next_available_aroma = 0.0;
        double projected_vol_loss = next_available_aroma * volatility_rate * dt_min;
        if (projected_vol_loss > next_available_aroma) projected_vol_loss = next_available_aroma;
        double next_vol_rate = projected_vol_loss / dt_min;
        if (next_vol_rate < 0.0) next_vol_rate = 0.0;

        projected_aroma_ext_rate = next_ext_rate;
        projected_aroma_vol_rate = next_vol_rate;
        projected_guard_margin = next_ext_rate - next_vol_rate;
        if (next_ext_rate > 1e-6) {
            projected_aroma_balance = next_vol_rate / (next_ext_rate + 1e-6);
        } else {
            projected_aroma_balance = (next_vol_rate > 0.0) ? 5.0 : 0.0;
        }
    }
    _Bool projected_flux_gate = (projected_aroma_ext_rate >= min_flux_ext && projected_aroma_vol_rate >= min_flux_vol);
    double aroma_balance_quality = 0.0;
    _Bool aroma_balance_cut = 0;
    double ratio_guard_ms = compute_ratio_diffusion_guard_ms(state);
    _Bool ratio_guard_met = (elapsed_ms >= ratio_guard_ms);
    if (!is_active_boiling) {
        double tolerance = get_aroma_balance_tolerance(state->tea_type);
        double percept_floor = aroma_perception_floor_ratio;
        double min_rate = get_aroma_signal_floor(state->tea_type);
        double min_elapsed_ms = get_aroma_balance_min_time_ms(state->tea_type);
        if (min_elapsed_ms < ratio_guard_ms) {
            min_elapsed_ms = ratio_guard_ms;
        }
        const double projection_threshold = 0.66;
        double balance_perception_gate = clampd(percept_floor + tolerance * 0.5, percept_floor, 0.9);
        if (aroma_capture_ratio >= balance_perception_gate &&
            live_flux_gate &&
            projected_aroma_ext_rate >= min_rate &&
            projected_aroma_vol_rate >= min_rate * 0.5 &&
            projected_flux_gate &&
            projected_aroma_balance >= projection_threshold &&
            projected_guard_margin <= 0.0 &&
            elapsed_ms >= min_elapsed_ms) {
            aroma_balance_cut = 1;
            double overshoot = projected_aroma_balance - projection_threshold;
            if (overshoot < 0.0) overshoot = 0.0;
            double normalized = overshoot / (tolerance + 1e-6);
            if (normalized > 1.0) normalized = 1.0;
            aroma_balance_quality = 1.0 - normalized;
        }
    }

    double amino_ratio = 0.0;
    if (cp.s_amino_acid > 0.0) {
        amino_ratio = clampd((e_ami_curr + de_ami) / cp.s_amino_acid, 0.0, 1.2);
    }
    set_br_d(&state->amino_depth_axis, amino_ratio);
    double amino_rate = (dt_min > 0.0) ? de_ami / dt_min : 0.0;
    set_br_d(&state->amino_vibrancy_axis, amino_rate);

    // 4. Analysis & Termination Logic
    double r_cat = (e_cat + de_cat) / cp.s_catechin;
    double r_caf = (e_caf + de_caf) / cp.s_caffeine;
    double r_ami = (e_ami + de_ami) / cp.s_amino_acid;
    double r_pol = (e_pol + de_pol) / cp.s_polysaccharide;
    double r_pec = (e_pec + de_pec) / cp.s_pectin;

    double sat_idx = (r_cat + r_ami + r_caf + r_pec + r_pol) / 5.0;
    set_br_d(&state->saturation_index, sat_idx);
    if (sat_idx > 0.99) state->is_exhausted = 1;

    double db_dt = (de_cat/cp.s_catechin * 0.7 + de_caf/cp.s_caffeine * 0.3) / dt_min;
    double ds_dt = (de_ami/cp.s_amino_acid * 0.4 + de_pol/cp.s_polysaccharide * 0.4 + de_pec/cp.s_pectin * 0.2) / dt_min;

    double prev_db_dt = get_double_br(&state->astringency_rate);
    
    double d2b_dt2 = (db_dt - prev_db_dt) / dt_min;
    double d2s_dt2 = (ds_dt - get_double_br(&state->sweetness_rate)) / dt_min;

    set_br_d(&state->astringency_rate, db_dt);
    set_br_d(&state->sweetness_rate, ds_dt);

    double override_progress = -1.0;
    _Bool style_override = 0;
    if (!is_active_boiling) {
        // Absolute aroma capture gates: allow an early stop once 0.03 mg is retained,
        // or once 0.02 mg is present and volatility rises to a 6:4 split.
        const double aroma_mass_full = 0.03;
        const double aroma_mass_min = 0.02;
        const double aroma_ratio_gate = 2.0 / 3.0; // 6:4 extraction-to-volatility
        double current_aroma_mass = aroma_retained;
        _Bool flux_present = (aroma_ext_rate > 1e-6 && aroma_vol_rate > 1e-6);
        _Bool aroma_mass_cut = 0;
        double aroma_mass_progress = -1.0;
        if (ratio_guard_met && flux_present && current_aroma_mass >= aroma_mass_full) {
            aroma_mass_cut = 1;
        } else if (ratio_guard_met && flux_present && current_aroma_mass >= aroma_mass_min &&
                   aroma_balance_ratio >= aroma_ratio_gate) {
            aroma_mass_cut = 1;
        }
        if (aroma_mass_cut) {
            double normalized_mass = clampd(current_aroma_mass / aroma_mass_full, 0.0, 1.25);
            aroma_mass_progress = normalized_mass * 100.0;
        }

        double aroma_safety_ratio = 0.0;
        if (aroma_ext_rate > 0.0) {
            aroma_safety_ratio = clampd(aroma_guard_margin / (aroma_ext_rate + 1e-6), -1.0, 1.0);
        }
        double aromatic_window_score = aroma_capture_ratio * 0.7 + (0.5 + 0.5 * aroma_safety_ratio) * 0.3;
        _Bool aromatic_window = (live_flux_gate && aroma_capture_ratio >= 0.62 && aroma_guard_margin > 0.003 && aromatic_window_score > 0.55);
        _Bool amino_window = (amino_ratio >= 0.35 && amino_rate >= 0.015);
        double balance_progress = (aroma_balance_cut)
            ? (aroma_capture_ratio * 0.6 + aroma_balance_quality * 0.4) * 100.0
            : -1.0;
        if (aroma_balance_cut || aroma_mass_cut) {
            style_override = 1;
            state->stop_signal = 1;
            double progress_floor = sat_idx * 100.0;
            double candidate_progress = -1.0;
            if (balance_progress >= 0.0) {
                candidate_progress = balance_progress;
            }
            if (aroma_mass_cut && aroma_mass_progress >= 0.0 &&
                (candidate_progress < 0.0 || aroma_mass_progress < candidate_progress)) {
                candidate_progress = aroma_mass_progress;
            }
            if (candidate_progress < progress_floor) candidate_progress = progress_floor;
            override_progress = candidate_progress;
        } else if (aromatic_window && amino_window && ratio_guard_met) {
            style_override = 1;
            state->stop_signal = 1;
            double aromatic_progress = (aroma_capture_ratio * 0.65 + amino_ratio * 0.35) * 100.0;
            double progress_floor = sat_idx * 100.0;
            if (aromatic_progress < progress_floor) aromatic_progress = progress_floor;
            override_progress = aromatic_progress;
        }

        double clarity_target = 1.3;
        if (heavy_solutes_conc > clarity_target) {
            style_override = 1;
            state->stop_signal = 1;
            double clarity_ratio = clampd(clarity_target / heavy_solutes_conc, 0.45, 1.0);
            double clarity_progress = sat_idx * 100.0 * clarity_ratio;
            if (override_progress < 0.0 || clarity_progress < override_progress) {
                override_progress = clarity_progress;
            }
        }
    }

    if (is_active_boiling) {
        // Decoction Termination: Viscosity / Total Dissolved Solids Concentration
        // Viscosity Index = (Pectin*1.5 + Polysaccharides*1.0 + Lignin*2.0) * Leaf Mass / Current Volume
        // This effectively models the "Thickening" as water evaporates.
        
        double decoction_conc = (e_pec * 1.5 + e_pol * 1.0 + e_lig * 2.0) * l_mass / current_vol;
        
        // Target roughly corresponds to "Thick" soup texture (e.g. > 3.0 mg/ml of heavy colloids)
        double target_viscosity = 3.0; 
        
        if (decoction_conc >= target_viscosity && elapsed_ms > 5000) {
            state->stop_signal = 1;
        } else {
            state->stop_signal = 0;
        }
        // Force extraction progress to follow viscosity metric in this mode but don't jump backwards
        double visc_progress = (decoction_conc / target_viscosity) * 100.0;
        double current_ext = get_double_br(&state->extraction_level);
        if (visc_progress < current_ext) visc_progress = current_ext;
        set_br_d(&state->extraction_level, visc_progress);
        
    } else if (!style_override) {
        // Gongfu / Standard Termination
        if (aroma_floor_met && ratio_guard_met && d2b_dt2 > d2s_dt2) {
            state->stop_signal = 1;
        } else {
            state->stop_signal = 0;
        }
        
        // Fix: Use saturation index for progress but clamp if stopped by signal
    }

    if (!is_active_boiling) {
        double progress = sat_idx * 100.0;
        if (style_override) {
            double final_progress = override_progress;
            if (final_progress < 0.0) final_progress = progress;
            set_br_d(&state->extraction_level, final_progress);
        } else {
            if (state->stop_signal) {
                // If stopped by bitterness, we might not have reached 100%
                // but we don't want to hang.
            }
            set_br_d(&state->extraction_level, progress);
        }
    }

    set_br_d(&state->ext_velocity, ds_dt * 100.0); 
    tea_distributor_refresh_snapshot(state, &tp, &cp);

    state->last_update_time = now;
}

double tea_distributor_adjust_target(tea_state_t *state, double base_slice_pct) {
    if (!state || base_slice_pct <= 0.0 || state->num_infusions == 0) {
        return base_slice_pct;
    }

    tea_profile_t tp;
    get_tea_profile(state->tea_type, &tp);
    tea_comp_profile_t cp;
    get_tea_comp_profile(state->tea_type, &cp);

    tea_distributor_refresh_snapshot(state, &tp, &cp);

    int infusions_left = (int)(state->num_infusions - state->current_infusion + 1);
    if (infusions_left < 1) infusions_left = 1;

    double remaining_pct = fmax(0.0, 100.0 - get_double_br(&state->extraction_level));
    double base_remaining = remaining_pct / infusions_left;
    if (base_slice_pct > base_remaining) base_slice_pct = base_remaining;

    double potential = get_double_br(&state->distributor_potential);
    if (potential <= 0.0) potential = 10.0;
    double conc_ratio = clampd(potential / 15.0, 0.5, 1.6);

    double length_factor = 1.0;
    compute_leaf_geometry_factor(state, &tp, &length_factor, NULL);
    double velocity_factor = clampd(tp.base_velocity / 0.12, 0.6, 1.5);

    double drag_ratio = get_double_br(&state->distributor_drag_ratio);
    if (drag_ratio <= 0.0) drag_ratio = 1.0;

    double distribution_factor = (conc_ratio + velocity_factor + length_factor) / 3.0;
    distribution_factor = clampd(distribution_factor * drag_ratio, 0.4, 1.6);

    double adjusted = base_slice_pct * distribution_factor;
    double min_slice = base_remaining * 0.35;
    if (adjusted < min_slice) adjusted = min_slice;
    if (adjusted > remaining_pct) adjusted = remaining_pct;

    set_br_d(&state->distributor_target_hint, adjusted);
    return adjusted;
}

// New helper to simulate lag
void physics_simulate_lag(tea_state_t *state, uint64_t lag_ms) {
    if (lag_ms == 0) return;
    
    // Step size 1 sec
    uint64_t step = 1000;
    uint64_t elapsed = 0;
    uint64_t start_vtime = state->last_update_time;
    
    while(elapsed < lag_ms) {
        uint64_t next_step = (lag_ms - elapsed > step) ? step : (lag_ms - elapsed);
        physics_simulate_step(state, start_vtime + elapsed + next_step);
        elapsed += next_step;
    }
}

void physics_get_extraction_rate(ttak_bigreal_t *rate_out, const tea_state_t *state) {
    ttak_bigreal_copy(rate_out, &state->ext_velocity, 0);
}

void physics_get_bitterness_accel(ttak_bigreal_t *accel_out, const tea_state_t *state) {
    ttak_bigreal_copy(accel_out, &state->bitter_accel, 0);
}

void tea_state_init(tea_state_t *state) {
    if (!state) return;
    
    uint64_t now = ttak_get_tick_count();
    
    #define INIT_BR(field) ttak_bigreal_init(&state->field, now);

    INIT_BR(current_temp)
    INIT_BR(extraction_level)
    INIT_BR(ext_catechin)
    INIT_BR(ext_amino_acid)
    INIT_BR(ext_caffeine)
    INIT_BR(ext_pectin)
    INIT_BR(ext_polysaccharide)
    INIT_BR(ext_lignin)
    INIT_BR(ext_aroma)
    INIT_BR(accum_ext_catechin)
    INIT_BR(accum_ext_amino_acid)
    INIT_BR(accum_ext_caffeine)
    INIT_BR(accum_ext_pectin)
    INIT_BR(accum_ext_polysaccharide)
    INIT_BR(accum_ext_lignin)
    INIT_BR(accum_ext_aroma)
    INIT_BR(structural_integrity)
    INIT_BR(hydration_state)
    INIT_BR(unfurling_state)
    INIT_BR(leaf_density)
    INIT_BR(astringency_rate)
    INIT_BR(sweetness_rate)
    INIT_BR(saturation_index)
    INIT_BR(ext_velocity)
    INIT_BR(bitter_accel)
    INIT_BR(aroma_extraction_axis)
    INIT_BR(aroma_volatility_axis)
    INIT_BR(amino_depth_axis)
    INIT_BR(amino_vibrancy_axis)
    INIT_BR(clarity_index)
    INIT_BR(distributor_potential)
    INIT_BR(distributor_drag_ratio)
    INIT_BR(distributor_target_hint)
    INIT_BR(water_hardness)
    INIT_BR(leaf_width)
    INIT_BR(leaf_height)
    INIT_BR(leaf_mass)
    INIT_BR(water_volume_ml)
    INIT_BR(altitude_m)
    INIT_BR(tds)
    INIT_BR(target_extraction_for_cycle)
    INIT_BR(cycle_start_extraction)

    for (int i = 0; i < 5; i++) {
        ttak_bigreal_init(&state->mineral_content[i], now);
    }
    #undef INIT_BR

    set_br_d(&state->current_temp, 0.0);
    set_br_d(&state->extraction_level, 0.0);
    set_br_d(&state->ext_catechin, 0.0);
    set_br_d(&state->ext_amino_acid, 0.0);
    set_br_d(&state->ext_caffeine, 0.0);
    set_br_d(&state->ext_pectin, 0.0);
    set_br_d(&state->ext_polysaccharide, 0.0);
    set_br_d(&state->ext_lignin, 0.0);
    set_br_d(&state->ext_aroma, 0.0);
    set_br_d(&state->accum_ext_catechin, 0.0);
    set_br_d(&state->accum_ext_amino_acid, 0.0);
    set_br_d(&state->accum_ext_caffeine, 0.0);
    set_br_d(&state->accum_ext_pectin, 0.0);
    set_br_d(&state->accum_ext_polysaccharide, 0.0);
    set_br_d(&state->accum_ext_lignin, 0.0);
    set_br_d(&state->accum_ext_aroma, 0.0);
    set_br_d(&state->structural_integrity, 1.0);
    set_br_d(&state->hydration_state, 0.0);
    set_br_d(&state->unfurling_state, 0.0);
    set_br_d(&state->leaf_density, 0.0);
    set_br_d(&state->astringency_rate, 0.0);
    set_br_d(&state->sweetness_rate, 0.0);
    set_br_d(&state->saturation_index, 0.0);
    set_br_d(&state->ext_velocity, 0.0);
    set_br_d(&state->bitter_accel, 0.0);
    set_br_d(&state->aroma_extraction_axis, 0.0);
    set_br_d(&state->aroma_volatility_axis, 0.0);
    set_br_d(&state->amino_depth_axis, 0.0);
    set_br_d(&state->amino_vibrancy_axis, 0.0);
    set_br_d(&state->clarity_index, 0.0);
    set_br_d(&state->distributor_potential, 0.0);
    set_br_d(&state->distributor_drag_ratio, 1.0);
    set_br_d(&state->distributor_target_hint, 0.0);
    set_br_d(&state->water_hardness, 0.0);
    set_br_d(&state->leaf_width, 0.0);
    set_br_d(&state->leaf_height, 0.0);
    set_br_d(&state->leaf_mass, 0.0);
    set_br_d(&state->water_volume_ml, 0.0);
    set_br_d(&state->altitude_m, 0.0);
    set_br_d(&state->tds, 0.0);
    set_br_d(&state->target_extraction_for_cycle, 0.0);
    set_br_d(&state->cycle_start_extraction, 0.0);
    state->temp_preserved = 0;
    state->is_exhausted = 0;
    for (int i = 0; i < 5; i++) {
        set_br_d(&state->mineral_content[i], 0.0);
    }
}

void tea_state_cleanup(void *ptr) {
    if (!ptr) return;
    tea_state_t *state = (tea_state_t *)ptr;
    ttak_bigreal_free(&state->current_temp, 0);
    ttak_bigreal_free(&state->extraction_level, 0);
    ttak_bigreal_free(&state->ext_catechin, 0);
    ttak_bigreal_free(&state->ext_amino_acid, 0);
    ttak_bigreal_free(&state->ext_caffeine, 0);
    ttak_bigreal_free(&state->ext_pectin, 0);
    ttak_bigreal_free(&state->ext_polysaccharide, 0);
    ttak_bigreal_free(&state->ext_lignin, 0);
    ttak_bigreal_free(&state->ext_aroma, 0);
    ttak_bigreal_free(&state->accum_ext_catechin, 0);
    ttak_bigreal_free(&state->accum_ext_amino_acid, 0);
    ttak_bigreal_free(&state->accum_ext_caffeine, 0);
    ttak_bigreal_free(&state->accum_ext_pectin, 0);
    ttak_bigreal_free(&state->accum_ext_polysaccharide, 0);
    ttak_bigreal_free(&state->accum_ext_lignin, 0);
    ttak_bigreal_free(&state->accum_ext_aroma, 0);
    ttak_bigreal_free(&state->structural_integrity, 0);
    ttak_bigreal_free(&state->hydration_state, 0);
    ttak_bigreal_free(&state->unfurling_state, 0);
    ttak_bigreal_free(&state->leaf_density, 0);
    ttak_bigreal_free(&state->astringency_rate, 0);
    ttak_bigreal_free(&state->sweetness_rate, 0);
    ttak_bigreal_free(&state->saturation_index, 0);
    ttak_bigreal_free(&state->ext_velocity, 0);
    ttak_bigreal_free(&state->bitter_accel, 0);
    ttak_bigreal_free(&state->aroma_extraction_axis, 0);
    ttak_bigreal_free(&state->aroma_volatility_axis, 0);
    ttak_bigreal_free(&state->amino_depth_axis, 0);
    ttak_bigreal_free(&state->amino_vibrancy_axis, 0);
    ttak_bigreal_free(&state->clarity_index, 0);
    ttak_bigreal_free(&state->distributor_potential, 0);
    ttak_bigreal_free(&state->distributor_drag_ratio, 0);
    ttak_bigreal_free(&state->distributor_target_hint, 0);
    ttak_bigreal_free(&state->water_hardness, 0);
    ttak_bigreal_free(&state->leaf_width, 0);
    ttak_bigreal_free(&state->leaf_height, 0);
    ttak_bigreal_free(&state->leaf_mass, 0);
    ttak_bigreal_free(&state->water_volume_ml, 0);
    ttak_bigreal_free(&state->altitude_m, 0);
    ttak_bigreal_free(&state->tds, 0);
    ttak_bigreal_free(&state->target_extraction_for_cycle, 0);
    ttak_bigreal_free(&state->cycle_start_extraction, 0);
    for (int i = 0; i < 5; i++) {
        ttak_bigreal_free(&state->mineral_content[i], 0);
    }
    // Note: self is freed by ttak_shared_destroy
}

void physics_backtrack_cooling(ttak_bigreal_t *temp_out, ttak_bigreal_t *start_temp, uint64_t duration_ms) {
    // Simple exponential decay: T = T0 * exp(-k*t)
    // k ~ 0.02 per min
    double t_min = duration_ms / 60000.0;
    double t0 = get_double_br(start_temp);
    double t_final = t0 * exp(-0.02 * t_min);
    set_br_d(temp_out, t_final);
}

void physics_backtrack_extraction(ttak_bigreal_t *extract_out, const ttak_bigreal_t *temp_profile, uint64_t duration_ms) {
    // This is now handled by physics_simulate_lag in the TUI
    (void)temp_profile;
    (void)extract_out;
    (void)duration_ms;
}
