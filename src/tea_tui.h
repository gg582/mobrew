#ifndef TEA_TUI_H
#define TEA_TUI_H

#include "tea_state.h"

void tui_init(void);
void tui_cleanup(void);

// Returns 0 on success, -1 if cancelled
int tui_show_config_wizard(ttak_shared_tea_state_t *shared_state, ttak_owner_t *owner);

int ask_choice(const char *question, const char **options, int count);

void tui_run_dashboard(ttak_shared_tea_state_t *shared_state, ttak_owner_t *owner);

#endif
