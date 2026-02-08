#include "tea_tui.h"
#include "tea_physics.h"
#include <ncurses.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <locale.h>
#include <ttak/timing/timing.h>
#include <ttak/stats/stats.h>
#include <ttak/mem/mem.h>

void tui_init(void) {
    setlocale(LC_ALL, "");
    initscr();
    cbreak();
    noecho();
    keypad(stdscr, TRUE);
    start_color();
    use_default_colors();
    init_pair(1, COLOR_GREEN, -1);
    init_pair(2, COLOR_CYAN, -1);
    init_pair(3, COLOR_RED, -1);
    init_pair(4, COLOR_YELLOW, -1); // For Golden Time
}

void tui_cleanup(void) {
    endwin();
}

static void print_centered(int row, const char *msg) {
    int y, x;
    getmaxyx(stdscr, y, x);
    (void)y;
    int start_x = (x - (int)strlen(msg)) / 2;
    if (start_x < 0) start_x = 0;
    mvprintw(row, start_x, "%s", msg);
}

// Helper to ask a question and get a choice
int ask_choice(const char *question, const char **options, int count) {
    int choice = 0;
    while(1) {
        clear();
        print_centered(2, question);
        for(int i=0; i<count; i++) {
            if(i == choice) attron(A_REVERSE);
            mvprintw(5+i, 4, "[%d] %s", i+1, options[i]);
            if(i == choice) attroff(A_REVERSE);
        }
        refresh();
        int ch = getch();
        if(ch == KEY_UP && choice > 0) choice--;
        if(ch == KEY_DOWN && choice < count-1) choice++;
        if(ch == 10) return choice; // Enter
    }
}

// Helper to get number input
static double ask_number(const char *question) {
    echo();
    char buf[64];
    if (question && strlen(question) > 0) {
        clear();
        print_centered(2, question);
        move(4, 4);
    } else {
        // Move to next line for clean micro-prompt
        int y, x;
        getyx(stdscr, y, x);
        (void)x; // Avoid unused variable warning
        move(y + 1, 4);
    }
    refresh();
    getnstr(buf, sizeof(buf) - 1);
    noecho();
    return (strlen(buf) > 0) ? atof(buf) : 0.0;
}

int tui_show_config_wizard(ttak_shared_tea_state_t *shared_state, ttak_owner_t *owner) {
    ttak_shared_result_t res;
    tea_state_t *state = ttak_shared_tea_state_access(shared_state, owner, &res);
    
    if(!state) {
        // Log the reason for failure if possible
        return -1;
    }
    
    // 0. Unit Preferences
    const char *v_units[] = {"Metric (ml, g)", "US Customary (US fl oz, oz)", "UK Imperial (UK fl oz, oz)"};
    state->unit_system = (unit_system_t)ask_choice("Select Volume/Mass Units:", v_units, 3);

    const char *t_units[] = {"Celsius (C)", "Fahrenheit (F)"};
    state->use_fahrenheit = (ask_choice("Select Temperature Units:", t_units, 2) == 1);

    const char *h_units[] = {"TDS (ppm, mg/L)", "French Degrees (fH)", "German Degrees (dH)"};
    state->hardness_unit = (hardness_unit_t)ask_choice("Select Water Hardness Units:", h_units, 3);

    // 0. Tea Type (New)
    const char *teas[] = {
        "Green (Normal)", "White Tea", "Black Tea", "Oolong Tea", "Yellow Tea",
        "Pu-Erh Tea", "Green (Gyokuro)", "Green (Sencha)", "Green (Fukamushi)", "Tibetan Fermented"
    };
    state->tea_type = (tea_type_t)ask_choice("Select Tea Profile:", teas, 10);

    // Tibetan specific
    state->has_butter = false;
    state->boiling_in_pot = false;
    state->heat_level = 0;
    state->is_ripe = false;

    if (state->tea_type == TEA_PUERH) {
        const char *pu_types[] = {"Sheng (Raw/Green)", "Shu (Ripe/Black)"};
        state->is_ripe = (ask_choice("Select Pu-Erh Type:", pu_types, 2) == 1);
    }

    // Vintage
    double years = ask_number("Vintage (Years since production):");
    state->vintage_years = (uint32_t)years;

    if (state->tea_type == TEA_TIBETAN || (state->tea_type == TEA_PUERH && state->vintage_years >= 15)) {
        const char *yn[] = {"Yes", "No"};
        bool boil = false;
        
        if (state->tea_type == TEA_TIBETAN) {
            if (ask_choice("Are you adding butter for thermal insulation? (Y/N)", yn, 2) == 0) {
                state->has_butter = true;
                boil = true;
            } else {
                if (ask_choice("Are you boiling the tea in a pot? (Y/N)", yn, 2) == 0) {
                    boil = true;
                }
            }
        } else {
            // Aged Pu-erh
            if (ask_choice("Aged Pu-Erh can be boiled. Are you boiling it in a pot? (Y/N)", yn, 2) == 0) {
                boil = true;
            }
        }
        
        if (boil) {
            state->boiling_in_pot = true;
            state->heat_level = (uint8_t)ask_number("How high did you turn the heat on in the pot? (1 is the weakest, 10 is the strongest)");
        }
    }

    // Leaf Geometry
    double w = ask_number("Leaf Width (Fully Wet) (mm):");
    double h = ask_number("Leaf Height (Fully Wet) (mm):");
    set_br_d(&state->leaf_width, w);
    set_br_d(&state->leaf_height, h);

    // Number of Infusions
    state->num_infusions = (uint32_t)ask_number("Number of Infusions:");
    if (state->num_infusions < 1) state->num_infusions = 1;
    state->current_infusion = 0; 
    state->cycle_active = false;

    // Leaf Mass & Water Volume
    char mass_prompt[64], vol_prompt[64];
    const char *v_unit_label = "ml";
    const char *m_unit_label = "g";
    
    if (state->unit_system == UNIT_US) { v_unit_label = "US fl oz"; m_unit_label = "oz"; }
    else if (state->unit_system == UNIT_UK) { v_unit_label = "UK fl oz"; m_unit_label = "oz"; }

    snprintf(mass_prompt, sizeof(mass_prompt), "Tea Mass (%s):", m_unit_label);
    snprintf(vol_prompt, sizeof(vol_prompt), "Water Volume (%s):", v_unit_label);
    
    double mass = ask_number(mass_prompt);
    double volume = ask_number(vol_prompt);
    
    if (state->unit_system != UNIT_METRIC) {
        mass *= 28.3495; // oz to g (common for both US/UK mass)
        if (state->unit_system == UNIT_US) {
            volume *= 29.5735; // US fl oz to ml
        } else {
            volume *= 28.4131; // UK fl oz to ml
        }
    }
    
    set_br_d(&state->leaf_mass, mass);
    set_br_d(&state->water_volume_ml, volume);

    // 1. Vessel
    const char *vessels[] = {
        "British Teapot", 
        "Zisha (Purple Clay)", 
        "White Porcelain (Teapot)", 
        "Celadon (Teapot)", 
        "Japanese Ceramics", 
        "Glass", 
        "White Porcelain Gaiwan", 
        "Celadon Gaiwan"
    };
    state->vessel = (vessel_type_t)ask_choice("Select Vessel Material:", vessels, 8);

    // 2. Boiling Method
    const char *boils[] = {
        "Electric Kettle",
        "Traditional Pot (Iron)",
        "Traditional Pot (Bronze)",
        "Traditional Pot (Clay)"
    };
    state->boil_method = (boil_method_t)ask_choice("Select Boiling Method:", boils, 4);

    // 2.5 Water Quality Wizard
    const char *yn_choices[] = {"Yes", "No"};
    state->is_bottled = (ask_choice("Are you using Bottled Mineral Water? (Y/N)", yn_choices, 2) == 0);
    
    if (state->is_bottled) {
        clear();
        print_centered(2, "=== Bottled Water Mineral Composition ===");
        mvprintw(5, 4, "Please enter the mineral content (mg/L) from the label:");
        double ca = ask_number("Calcium (Ca):");
        double mg = ask_number("Magnesium (Mg):");
        double na = ask_number("Sodium (Na):");
        double k  = ask_number("Potassium (K):");
        set_br_d(&state->mineral_content[0], ca);
        set_br_d(&state->mineral_content[1], mg);
        set_br_d(&state->mineral_content[2], na);
        set_br_d(&state->mineral_content[3], k);
        // Calculate approx TDS from minerals
        set_br_d(&state->tds, ca + mg + na + k);
    } else {
        const char *know_tds[] = {"Yes, I know the hardness value", "No, I don't know"};
        if (ask_choice("Do you know the hardness of your tap water?", know_tds, 2) == 0) {
            char h_prompt[64];
            const char *h_label = "ppm";
            if (state->hardness_unit == HARDNESS_FH) h_label = "fH";
            else if (state->hardness_unit == HARDNESS_DH) h_label = "dH";
            
            snprintf(h_prompt, sizeof(h_prompt), "Enter hardness value (%s):", h_label);
            double h_val = ask_number(h_prompt);
            
            // Normalize to ppm
            if (state->hardness_unit == HARDNESS_FH) h_val *= 10.0;
            else if (state->hardness_unit == HARDNESS_DH) h_val *= 17.848;
            
            set_br_d(&state->tds, h_val);
        } else {
            // User doesn't know TDS
            const char *bedrock[] = {"Yes (Hard water likely)", "No (Soft water likely)"};
            state->is_limestone_bedrock = (ask_choice("Is the local bedrock in your area limestone? (Y/N)", bedrock, 2) == 0);
            // Assign default TDS based on bedrock
            set_br_d(&state->tds, state->is_limestone_bedrock ? 250.0 : 50.0);
        }
    }

    // Altitude
    double altitude = ask_number("If you live in a mountainous or lowland area, enter your altitude (m), otherwise press Enter:");
    set_br_d(&state->altitude_m, altitude);
    double boiling_point = 100.0 - (altitude / 285.0);
    if (boiling_point < 70.0) boiling_point = 70.0;

    // Temperature Input
    double restored_temp = get_double_br(&state->current_temp);
    _Bool skip_temp_input = 0;
    if (state->temp_preserved && restored_temp > 20.0) {
        char restored_msg[128];
        snprintf(restored_msg, sizeof(restored_msg), "Restored temperature: %.1f %s. Use it?",
                 state->use_fahrenheit ? (restored_temp * 9.0/5.0 + 32.0) : restored_temp,
                 state->use_fahrenheit ? "F" : "C");
        const char *yn[] = {"Yes, use restored", "No, enter new"};
        if (ask_choice(restored_msg, yn, 2) == 0) {
            skip_temp_input = 1;
        }
    }

    if (!skip_temp_input) {
        clear();
        char temp_prompt[128];
        double disp_bp = state->use_fahrenheit ? (boiling_point * 9.0/5.0 + 32.0) : boiling_point;
        snprintf(temp_prompt, sizeof(temp_prompt), "Enter current water temperature (%s):", 
                 state->use_fahrenheit ? "F" : "C");
        
        print_centered(2, temp_prompt);
        mvprintw(3, 4, "(Altitude-adjusted Boiling Point: %.1f %s)", 
                 disp_bp, state->use_fahrenheit ? "F" : "C");
        
        double current_temp = ask_number("");
        if (state->use_fahrenheit) {
            current_temp = (current_temp - 32.0) * 5.0 / 9.0; // F to C
        }

        if (current_temp > boiling_point) {
            current_temp = boiling_point;
            mvprintw(6, 4, "Warning: Temperature capped at boiling point (%.1f C)", boiling_point);
            refresh();
            napms(1000);
        }
        set_br_d(&state->current_temp, current_temp);
    }

    // 3. Operational Lag
    const char *scenarios[] = {"Pre-poured (Water is in pot)", "In-Kettle (Just boiled / Boiling)"};
    int scenario = ask_choice("Current Status:", scenarios, 2);

    double lag_minutes = 0;

    if (scenario == 0) { // Pre-poured (Case A)
        lag_minutes = ask_number("How many minutes passed since pouring?");
    } else { // In-Kettle (Case B)
        const char *immediate[] = {"Yes, immediately after boil", "No, it boiled a while ago"};
        int imm = ask_choice("Did you enter this prompt immediately after boil?", immediate, 2);
        
        if (imm == 0) {
            lag_minutes = 0;
        } else {
            lag_minutes = ask_number("How many minutes passed since it boiled?");
        }
    }

    // Initialize Physics Model State
    set_br_d(&state->extraction_level, 0.0);
    set_br_d(&state->ext_catechin, 0.0);
    set_br_d(&state->ext_amino_acid, 0.0);
    set_br_d(&state->ext_caffeine, 0.0);
    set_br_d(&state->ext_pectin, 0.0);
    set_br_d(&state->ext_polysaccharide, 0.0);
    set_br_d(&state->ext_lignin, 0.0);
    set_br_d(&state->accum_ext_catechin, 0.0);
    set_br_d(&state->accum_ext_amino_acid, 0.0);
    set_br_d(&state->accum_ext_caffeine, 0.0);
    set_br_d(&state->accum_ext_pectin, 0.0);
    set_br_d(&state->accum_ext_polysaccharide, 0.0);
    set_br_d(&state->accum_ext_lignin, 0.0);
    set_br_d(&state->hydration_state, 0.0);
    set_br_d(&state->structural_integrity, 1.0); 
    set_br_d(&state->astringency_rate, 0.0);

    tea_profile_t tp;
    get_tea_profile(state->tea_type, &tp);
    set_br_d(&state->leaf_density, tp.leaf_density);
    
    state->start_time = ttak_get_tick_count();
    state->last_update_time = state->start_time;
    state->cycle_start_time = state->start_time;

    // Physics Backtracking / Simulation
    if (lag_minutes > 0) {
        uint64_t lag_ms = (uint64_t)(lag_minutes * 60 * 1000);
        
        if (scenario == 0) {
            // Case A: Water on leaves -> Extraction happened!
            // Step-by-step simulation for accuracy
             state->last_update_time = state->start_time - lag_ms;
             state->cycle_start_time = state->last_update_time;
             physics_simulate_lag(state, lag_ms); 
        } else {
            // Case B: Water in kettle -> Just Cooling
            ttak_bigreal_t temp_out;
            ttak_bigreal_init(&temp_out, 0);
            physics_backtrack_cooling(&temp_out, &state->current_temp, lag_ms);
            ttak_bigreal_copy(&state->current_temp, &temp_out, 0);
            ttak_bigreal_free(&temp_out, 0);
        }
    }

    ttak_shared_tea_state_release(shared_state);
    return 0;
}

static void show_consistency_stats(void) {
    clear();
    print_centered(2, "=== Extraction Consistency (Simulated History) ===");
    
    // Generate dummy data (extraction percentages * 100)
    size_t count = 100;
    uint64_t now = ttak_get_tick_count();
    uint64_t *data = ttak_mem_alloc(count * sizeof(uint64_t), 10000, now);
    if (!data) return;
    
    for(size_t i=0; i<count; i++) {
        data[i] = 80 + (rand() % 40); // 80 - 120%
    }
    
    ttak_bigreal_t p50, p95, p99, p999;
    
    ttak_bigreal_init(&p50, now);
    ttak_bigreal_init(&p95, now);
    ttak_bigreal_init(&p99, now);
    ttak_bigreal_init(&p999, now);
    
    ttak_stats_compute_percentiles(data, count, &p50, &p95, &p99, &p999, now);
    
    uint64_t val50=0, val95=0, val99=0;
    ttak_bigint_export_u64(&p50.mantissa, &val50);
    ttak_bigint_export_u64(&p95.mantissa, &val95);
    ttak_bigint_export_u64(&p99.mantissa, &val99);
    
    mvprintw(5, 5, "P50 (Median): %llu %%", (unsigned long long)val50);
    mvprintw(6, 5, "P95         : %llu %%", (unsigned long long)val95);
    mvprintw(7, 5, "P99         : %llu %%", (unsigned long long)val99);
    
    mvprintw(10, 5, "Press any key to return...");
    refresh();
    
    nodelay(stdscr, FALSE);
    getch();
    nodelay(stdscr, TRUE);
}

typedef struct {
    float x, y;
    float vx, vy;
    bool active;
} particle_t;

#define MAX_PARTICLES 200
static particle_t particles[MAX_PARTICLES];
static bool particles_init = false;

static void update_particles(double velocity, double level) {
    if (!particles_init) {
        for (int i = 0; i < MAX_PARTICLES; i++) {
            particles[i].active = false;
        }
        particles_init = true;
    }

    // Number of active particles based on level (0-100)
    int target_active = (int)(level * 2.0);
    if (target_active > MAX_PARTICLES) target_active = MAX_PARTICLES;
    if (target_active < 5 && level > 0.1) target_active = 5;

    // Base speed from velocity (ds_dt * 100)
    float speed_factor = (float)velocity * 0.05f;
    if (speed_factor < 0.05f && level > 0) speed_factor = 0.05f;
    if (speed_factor > 2.0f) speed_factor = 2.0f;

    for (int i = 0; i < MAX_PARTICLES; i++) {
        if (i < target_active) {
            if (!particles[i].active) {
                particles[i].active = true;
                particles[i].x = (float)(rand() % 38 + 1);
                particles[i].y = (float)(rand() % 8 + 1);
                particles[i].vx = ((float)(rand() % 100) / 50.0f - 1.0f) * speed_factor;
                particles[i].vy = ((float)(rand() % 100) / 50.0f - 1.0f) * speed_factor;
            } else {
                particles[i].x += particles[i].vx;
                particles[i].y += particles[i].vy;

                // Decay velocity slightly for organic feel
                particles[i].vx *= 0.99f;
                particles[i].vy *= 0.99f;

                // Bounce
                if (particles[i].x <= 1 || particles[i].x >= 39) particles[i].vx *= -1;
                if (particles[i].y <= 1 || particles[i].y >= 9) particles[i].vy *= -1;
            }
        } else {
            particles[i].active = false;
        }
    }
}

static void draw_leaching_box(int start_y, int start_x, double level, double velocity) {
    // Unicode box using thick lines
    attron(COLOR_PAIR(2));
    mvprintw(start_y, start_x, "┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓");
    for (int i = 1; i < 10; i++) {
        mvprintw(start_y + i, start_x, "┃                                        ┃");
    }
    mvprintw(start_y + 10, start_x, "┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛");
    attroff(COLOR_PAIR(2));

    // Choose shading block based on velocity and concentration
    const char *shade = "░";
    if (velocity > 0.6) shade = "⣿";
    else if (velocity > 0.3) shade = "▓";
    else if (level > 60.0) shade = "▒";
    else if (level > 20.0) shade = "░";
    else if (level > 0.1) shade = "·"; // Early diffusion

    attron(COLOR_PAIR(1));
    for (int i = 0; i < MAX_PARTICLES; i++) {
        if (particles[i].active) {
            // Add some jitter to the shade for a more organic fluid look
            const char *p_shade = shade;
            if (velocity > 0.1 && (rand() % 5 == 0)) p_shade = "▒";
            
            mvprintw(start_y + (int)particles[i].y, start_x + (int)particles[i].x, "%s", p_shade);
        }
    }
    attroff(COLOR_PAIR(1));
    
    mvprintw(start_y, start_x + 2, "[ Tea Leaching Simulator ]");
}

static void draw_leaf_render(int start_y, int start_x, double unfurl) {
    attron(COLOR_PAIR(1) | A_BOLD);
    mvprintw(start_y, start_x, "Leaf Morphology (Unicode Block Art)");
    
    if (unfurl < 0.25) {
        // Stage 1: Tightly Rolled / Pearl style
        mvprintw(start_y + 2, start_x + 5, "   ⣴⣾⣦   ");
        mvprintw(start_y + 3, start_x + 5, "   ⢿⣿⡿   ");
        mvprintw(start_y + 4, start_x + 5, "   ⠻⢿⠟   ");
    } else if (unfurl < 0.50) {
        // Stage 2: Initial Swelling
        mvprintw(start_y + 2, start_x + 4, "  ⣠⣾⣿⣷⣄  ");
        mvprintw(start_y + 3, start_x + 4, " ⢀⣿⣿⣿⣿⡀ ");
        mvprintw(start_y + 4, start_x + 4, "  ⠙⢿⣿⡿⠋  ");
    } else if (unfurl < 0.75) {
        // Stage 3: Expanding Ellipse
        mvprintw(start_y + 1, start_x + 3, "   ⣠⣴⣾⣷⣦⣄   ");
        mvprintw(start_y + 2, start_x + 3, " ⣴⣿⣿⣿⣿⣿⣿⣦ ");
        mvprintw(start_y + 3, start_x + 3, " ⠻⣿⣿⣿⣿⣿⣿⠟ ");
        mvprintw(start_y + 4, start_x + 3, "   ⠙⠻⢿⡿⠟⠋   ");
    } else {
        // Stage 4: Fully Unfurled (Leaf with Vein and Stem)
        mvprintw(start_y + 1, start_x + 2, "      ⣴⣾⣦      ");
        mvprintw(start_y + 2, start_x + 2, "    ⣠⣾⣿⣿⣷⣄    ");
        mvprintw(start_y + 3, start_x + 2, "  ⣴⣿⣿⢿⣿⡿⣿⣿⣦  ");
        mvprintw(start_y + 4, start_x + 2, "   ⠙⢿⣿⣿⣿⣿⡿⠋   ");
        mvprintw(start_y + 5, start_x + 2, "      ⣿⣿      "); 
        mvprintw(start_y + 6, start_x + 2, "      ⠙⠋      ");
    }
    attroff(COLOR_PAIR(1) | A_BOLD);
}

void tui_run_dashboard(ttak_shared_tea_state_t *shared_state, ttak_owner_t *owner) {
    ttak_shared_result_t res;
    
    nodelay(stdscr, TRUE); 

    while(1) {
        uint64_t now = ttak_get_tick_count();
        tt_autoclean_dirty_pointers(now);
        
        tea_state_t *state = ttak_shared_tea_state_access(shared_state, owner, &res);
        if (state) {
            // Check if we need to start a new infusion cycle
            if (!state->cycle_active && state->current_infusion < state->num_infusions) {
                bool is_initial = (state->current_infusion == 0);
                
                // If one just finished, we wait for a moment or a key before jumping to next prompt
                // BUT if it's the very first one, we can start it.
                if (is_initial) {
                    state->current_infusion = 1;
                    
                    // Set extraction targets for the first cycle
                    ttak_bigreal_copy(&state->cycle_start_extraction, &state->extraction_level, 0);
                    set_br_d(&state->target_extraction_for_cycle, 100.0 / state->num_infusions);
                    
                    state->cycle_active = true;
                    // Reset cycle start time to now to ensure time starts from 0 for 1st infusion
                    state->cycle_start_time = ttak_get_tick_count();

                    // Ensure last_update_time is caught up to now
                    state->last_update_time = state->cycle_start_time;
                } else {
                    // This is a transition between infusions.
                    // We wait for the user to be ready so they can see the previous result.
                    
                    // Pause and ask for parameters
                    nodelay(stdscr, FALSE);
                    clear();
                    char cycle_msg[128];
                    snprintf(cycle_msg, sizeof(cycle_msg), "=== Preparing Infusion %u of %u ===", 
                             state->current_infusion + 1, state->num_infusions);
                    print_centered(2, cycle_msg);
                    mvprintw(4, 4, "Previous infusion complete. Ready for next?");
                    mvprintw(6, 4, "Press any key to enter parameters...");
                    refresh();
                    getch();

                    // Before starting next infusion, move current extraction to accumulated
                    ttak_bigreal_t tmp, total;
                    ttak_bigreal_init(&tmp, 0);
                    ttak_bigreal_init(&total, 0);
                    
                    #define ACCUM_EXT(field) \
                        ttak_bigreal_copy(&tmp, &state->accum_ext_##field, 0); \
                        ttak_bigreal_add(&total, &tmp, &state->ext_##field, 0); \
                        ttak_bigreal_copy(&state->accum_ext_##field, &total, 0); \
                        set_br_d(&state->ext_##field, 0.0);

                    ACCUM_EXT(catechin)
                    ACCUM_EXT(amino_acid)
                    ACCUM_EXT(caffeine)
                    ACCUM_EXT(pectin)
                    ACCUM_EXT(polysaccharide)
                    ACCUM_EXT(lignin)
                    #undef ACCUM_EXT
                    
                    ttak_bigreal_free(&tmp, 0);
                    ttak_bigreal_free(&total, 0);

                    // Reset derivatives and signals for fresh start
                    set_br_d(&state->astringency_rate, 0.0);
                    set_br_d(&state->sweetness_rate, 0.0);
                    set_br_d(&state->ext_velocity, 0.0);
                    set_br_d(&state->bitter_accel, 0.0);
                    state->stop_signal = false;

                    state->current_infusion++;
                    clear();
                    print_centered(2, cycle_msg);
                    
                    // Prompt for Temperature
                    double boiling_point = 100.0 - (get_double_br(&state->altitude_m) / 285.0);
                    if (boiling_point < 70.0) boiling_point = 70.0;
                    
                    char temp_prompt[128];
                    snprintf(temp_prompt, sizeof(temp_prompt), "Enter water temperature for this infusion (%s):", 
                             state->use_fahrenheit ? "F" : "C");
                    mvprintw(4, 4, "%s", temp_prompt);
                    refresh();
                    
                    double next_temp = ask_number("");
                    if (state->use_fahrenheit) next_temp = (next_temp - 32.0) * 5.0 / 9.0;
                    if (next_temp > boiling_point) next_temp = boiling_point;
                    set_br_d(&state->current_temp, next_temp);
                    
                    // Prompt for Volume
                    char vol_prompt[128];
                    const char *v_unit = "ml";
                    if (state->unit_system == UNIT_US) v_unit = "US fl oz";
                    else if (state->unit_system == UNIT_UK) v_unit = "UK fl oz";
                    
                    snprintf(vol_prompt, sizeof(vol_prompt), "Enter water volume for this infusion (%s):", v_unit);
                    mvprintw(6, 4, "%s", vol_prompt);
                    refresh();
                    
                    double next_vol = ask_number("");
                    if (state->unit_system == UNIT_US) next_vol *= 29.5735;
                    else if (state->unit_system == UNIT_UK) next_vol *= 28.4131;
                    set_br_d(&state->water_volume_ml, next_vol);
                    
                    // Set cycle start time
                    state->cycle_start_time = ttak_get_tick_count();
                    state->last_update_time = state->cycle_start_time;

                    // Set extraction targets
                    ttak_bigreal_copy(&state->cycle_start_extraction, &state->extraction_level, 0);
                    set_br_d(&state->target_extraction_for_cycle, (double)state->current_infusion * 100.0 / state->num_infusions);
                    
                    state->cycle_active = true;
                    nodelay(stdscr, TRUE);
                }
            }

            if (state->cycle_active) {
                physics_simulate_step(state, now);
                
                // Stop At Acceleration Surge: Signal brew completion when bitterness overwhelms sweetness
                if (state->stop_signal) {
                    state->cycle_active = false;
                }
                
                // Fallback: Check if target reached (Traditional saturation goal)
                double current_ext = get_double_br(&state->extraction_level);
                double target_ext = get_double_br(&state->target_extraction_for_cycle);
                
                if (current_ext >= target_ext) {
                    state->cycle_active = false;
                    set_br_d(&state->extraction_level, target_ext); // Clamp to target
                }
            }
            
            // Update particles based on extraction speed and total level
            update_particles(get_double_br(&state->ext_velocity), get_double_br(&state->extraction_level));

            clear();
            attron(COLOR_PAIR(1));
            print_centered(1, "=== MOBREW ===");
            attroff(COLOR_PAIR(1));

            mvprintw(3, 5, "Infusion: %u / %u", state->current_infusion, state->num_infusions);
            
            // Stopwatch / Brew Time
            if (state->cycle_active || get_double_br(&state->extraction_level) > 0) {
                uint64_t elapsed_ms = now - state->cycle_start_time;
                if (!state->cycle_active && state->last_update_time > state->cycle_start_time) {
                    elapsed_ms = state->last_update_time - state->cycle_start_time;
                }
                int emins = (int)(elapsed_ms / 1000 / 60);
                int esecs = (int)((elapsed_ms / 1000) % 60);
                attron(COLOR_PAIR(2) | A_BOLD);
                mvprintw(3, 35, "Brew Time:  %02d:%02d", emins, esecs);
                attroff(COLOR_PAIR(2) | A_BOLD);
            }

            if (!state->cycle_active) {
                if (state->is_exhausted) {
                    attron(COLOR_PAIR(3) | A_BOLD);
                    print_centered(18, "Leaf essence exhausted. No more extraction possible.");
                    print_centered(19, "I hope you enjoyed your tea time.");
                    attroff(COLOR_PAIR(3) | A_BOLD);
                }

                if (state->current_infusion == state->num_infusions) {
                    attron(COLOR_PAIR(3) | A_BOLD);
                    mvprintw(3, 22, "[ SESSIONS COMPLETE ]");
                    attroff(COLOR_PAIR(3) | A_BOLD);
                } else if (state->current_infusion > 0) {
                    mvprintw(3, 22, "[ INFUSION COMPLETE ]");
                }
            }

            // Temperature & extraction
            mvprintw(5, 5, "Temperature: %.1f %s", 
                     state->use_fahrenheit ? (get_double_br(&state->current_temp) * 9.0/5.0 + 32.0) : get_double_br(&state->current_temp),
                     state->use_fahrenheit ? "F" : "C");
            
            attron(A_BOLD);
            mvprintw(6, 5, "Balance Target : %.2f  \xCE\x94%+.3f/min", 
                     get_double_br(&state->saturation_index),
                     get_double_br(&state->ext_velocity) / 100.0);
            attroff(A_BOLD);

            mvprintw(7, 5, "Extraction : %.2f %%", get_double_br(&state->extraction_level));
            
            // Draw progress bar
            int bar_width = 40;
            double progress = get_double_br(&state->extraction_level) / 100.0;
            if (progress > 1.0) progress = 1.0;
            int filled = (int)(progress * bar_width);
            
            mvprintw(8, 5, "[");
            for(int i=0; i<bar_width; i++) {
                if(i < filled) addch('=');
                else addch(' ');
            }
            addch(']');

            // Detailed composition
            mvprintw(10, 5, "--- Composition ---");
            mvprintw(11, 7, "Catechin: %.2f mg", get_double_br(&state->ext_catechin));
            mvprintw(12, 7, "Theanine: %.2f mg", get_double_br(&state->ext_amino_acid));
            mvprintw(13, 7, "Caffeine: %.2f mg", get_double_br(&state->ext_caffeine));
            
            // Hydration Bar
            double h = get_double_br(&state->hydration_state);
            int h_filled = (int)round(h * 12.0);
            if (h_filled > 12) h_filled = 12;
            if (h_filled < 0) h_filled = 0;
            int h_empty = 12 - h_filled;
            mvprintw(14, 5, "Hydration : ");
            for (int i = 0; i < h_filled; i++) addstr("\xE2\x96\x88");
            for (int i = 0; i < h_empty; i++) addstr("\xE2\x96\x91");
            printw(" %.2f", h);
            
            // Saturation Index
            mvprintw(11, 35, "Saturation: %.1f %%", get_double_br(&state->saturation_index) * 100.0);
            
            // Draw Leaching Simulator Box (Updated call)
            draw_leaching_box(5, 50, get_double_br(&state->extraction_level), get_double_br(&state->ext_velocity) / 100.0);

            // Draw Leaf Render
            draw_leaf_render(16, 50, get_double_br(&state->unfurling_state));

            // Warning if over-extracted
            if (get_double_br(&state->astringency_rate) > 0.5) {
                attron(COLOR_PAIR(3) | A_BLINK);
                mvprintw(15, 5, "!!! WARNING: Bitterness Acceleration High !!!");
                attroff(COLOR_PAIR(3) | A_BLINK);
            }

            mvprintw(21, 5, "Commands: [S] Stats  [Q] Quit");
            
            refresh();
            ttak_shared_tea_state_release(shared_state);
        }

        int ch = getch();
        if (ch == 'q' || ch == 'Q') break;
        if (ch == 's' || ch == 'S') show_consistency_stats();
        
        napms(50); // Faster refresh for smoother particles
    }
}