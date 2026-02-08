#ifndef TEA_PERSISTENCE_H
#define TEA_PERSISTENCE_H

#include "tea_state.h"

void tea_save_state(const tea_state_t *state);
_Bool tea_load_state(tea_state_t *state, uint64_t *saved_timestamp);

#endif
