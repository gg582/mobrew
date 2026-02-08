#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ttak/timing/timing.h>
#include "tea_state.h"
#include "tea_tui.h"
#include "tea_physics.h"
#include "tea_persistence.h"

// Simple JSON-ish loader for presets
static int load_preset(const char *filename, tea_state_t *state) {
    FILE *f = fopen(filename, "r");
    if (!f) return -1;

    // Very basic key-value parser for the prototype
    char line[256];
    while (fgets(line, sizeof(line), f)) {
        char key[64], val[128];
        if (sscanf(line, " \"%63[^\"]\": \"%127[^\"]\"", key, val) == 2 ||
            sscanf(line, " \"%63[^\"]\": %127[^, \n]", key, val) == 2) {
            
            if (strcmp(key, "tea_type") == 0) state->tea_type = (tea_type_t)atoi(val);
            else if (strcmp(key, "vessel") == 0) state->vessel = (vessel_type_t)atoi(val);
            else if (strcmp(key, "leaf_mass") == 0) set_br_d(&state->leaf_mass, atof(val));
            else if (strcmp(key, "water_volume_ml") == 0) set_br_d(&state->water_volume_ml, atof(val));
            else if (strcmp(key, "current_temp") == 0) set_br_d(&state->current_temp, atof(val));
            else if (strcmp(key, "num_infusions") == 0) state->num_infusions = atoi(val);
            else if (strcmp(key, "tds") == 0) set_br_d(&state->tds, atof(val));
            else if (strcmp(key, "altitude_m") == 0) set_br_d(&state->altitude_m, atof(val));
        }
    }
    
    fclose(f);
    
    // Set up initial timing for preset mode
    state->start_time = ttak_get_tick_count();
    state->last_update_time = state->start_time;
    state->cycle_start_time = state->start_time;
    
    return 0;
}

int main(int argc, char **argv) {
    // 1. Create a safe owner
    ttak_owner_t *owner = ttak_owner_create(TTAK_OWNER_SAFE_DEFAULT);
    if (!owner) {
        fprintf(stderr, "Failed to create owner.\n");
        return 1;
    }

    // Parse arguments
    bool use_preset = false;
    const char *preset_file = NULL;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--preset") == 0 && i + 1 < argc) {
            preset_file = argv[i + 1];
            use_preset = true;
            i++;
        }
    }

    // 2. Setup Shared State
    ttak_shared_tea_state_t shared_state;
    ttak_shared_tea_state_init(&shared_state);
    shared_state.base.cleanup = tea_state_cleanup;
    
    // 3. Add owner to shared state
    if (shared_state.base.add_owner(&shared_state.base, owner) != TTAK_OWNER_SUCCESS) {
        fprintf(stderr, "Failed to add owner to shared state.\n");
        ttak_owner_destroy(owner);
        return 1;
    }

    // 4. Allocate payload
    if (ttak_shared_tea_state_allocate(&shared_state, TTAK_SHARED_LEVEL_3) != TTAK_OWNER_SUCCESS) {
        fprintf(stderr, "Failed to allocate shared tea state.\n");
        ttak_owner_destroy(owner);
        return 1;
    }

    // 5. Initialize payload members
    ttak_shared_result_t res;
    tea_state_t *state = ttak_shared_tea_state_access(&shared_state, owner, &res);
    uint64_t saved_ts = 0;
    if (state) {
        tea_state_init(state);
        
        if (use_preset) {
            if (load_preset(preset_file, state) != 0) {
                fprintf(stderr, "Error: Could not load preset file %s\n", preset_file);
                ttak_shared_tea_state_release(&shared_state);
                ttak_shared_destroy(&shared_state.base);
                ttak_owner_destroy(owner);
                return 1;
            }
        } else {
            if (tea_load_state(state, &saved_ts)) {
                state->temp_preserved = 1;
            }
        }
        ttak_shared_tea_state_release(&shared_state);
    } else {
        fprintf(stderr, "Failed to access state for initialization.\n");
        ttak_shared_destroy(&shared_state.base);
        ttak_owner_destroy(owner);
        return 1;
    }
    
    // 6. Init TUI
    tui_init();

    if (!use_preset) {
        // 7. Run Wizard
        if (tui_show_config_wizard(&shared_state, owner) != 0) {
            tui_cleanup();
            ttak_shared_destroy(&shared_state.base);
            ttak_owner_destroy(owner);
            return 0;
        }
    }
    
    // 8. Run Dashboard
    tui_run_dashboard(&shared_state, owner);
    
    // 8.5 Save state on exit
    state = ttak_shared_tea_state_access(&shared_state, owner, &res);
    if (state) {
        tea_save_state(state);
        ttak_shared_tea_state_release(&shared_state);
    }

    // 9. Cleanup
    tui_cleanup();
    
    ttak_shared_destroy(&shared_state.base);
    ttak_owner_destroy(owner);
    
    return 0;
}