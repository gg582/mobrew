#include "tea_persistence.h"
#include "tea_physics.h"
#include <stdio.h>
#include <stdlib.h>
#include <ttak/timing/timing.h>

#define STATE_FILE ".mobrew_state"

void tea_save_state(const tea_state_t *state) {
    if (!state) return;
    FILE *fp = fopen(STATE_FILE, "w");
    if (!fp) return;

    double temp = get_double_br(&state->current_temp);
    uint64_t now = ttak_get_tick_count();

    fprintf(fp, "%llu %f\n", (unsigned long long)now, temp);
    fclose(fp);
}

_Bool tea_load_state(tea_state_t *state, uint64_t *saved_timestamp) {
    if (!state) return 0;
    FILE *fp = fopen(STATE_FILE, "r");
    if (!fp) return 0;

    unsigned long long ts;
    double temp;
    if (fscanf(fp, "%llu %lf", &ts, &temp) == 2) {
        *saved_timestamp = (uint64_t)ts;
        set_br_d(&state->current_temp, temp);
        fclose(fp);
        return 1;
    }

    fclose(fp);
    return 0;
}