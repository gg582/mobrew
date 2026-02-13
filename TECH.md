# MOBREW: Universal Hydration & Extraction Simulator

`mobrew` is a high-fidelity tea brewing simulator that models the multi-stage physical and chemical kinetics of leaf expansion and compound dissolution. Unlike basic timers, it employs a **Universal Hydration Model** and **Cumulative Component Tracking** to provide real-time feedback on extraction balance and astringency limits.

---

## 1. Core Physics & Simulation Logic

### A. Hydration-Driven Geometry Model
*   **Concept**: Dry tea leaves have minimal effective surface area. As they absorb water, they expand toward their asymptotic limit ($A_{wet}$), increasing the contact area for extraction.
*   **Hydration State ($H$)**: A dimensionless value from 0.0 (bone dry) to 1.0 (fully saturated).
*   **Effective Surface Area ($A_{eff}$)**: Modeled as $A_{eff} = A_{wet} \times \sqrt{H}$. This ensures that even partially hydrated leaves contribute significantly to extraction while maintaining the initial "wetting delay."

### B. Multi-Round Persistence
*   **Logic**: Subsequent infusions do not reset the leaf's physical state. 
*   **State Persistence**: Leaves inherit their physical state (hydration, structural integrity) from the previous round naturally.
*   **Concentration Gradient**: Component extraction is cumulative. The drive for the next round is $S_i - (e_{i,accum} + e_{i,curr})$, where $S_i$ is the total soluble mass. Current session extraction ($e_{i,curr}$) starts from zero for each infusion.

### C. Non-Linear Component Kinetics
*   **Arrhenius Scaling**: Extraction rates ($k_i$) for individual components (Catechins, Amino Acids, Caffeine, etc.) scale exponentially with temperature.
*   **Sensitivity Factors**: Different components have unique temperature sensitivities. Amino acids (Umami) extract even at lower temperatures, while Catechins (Bitterness) require high kinetic energy.

---

## 2. Mathematical Models

### 2.1. Leaf Hydration Model
The hydration state $H$ evolves as a function of temperature ($T$), leaf density ($\rho$), and current saturation:
$$\frac{dH}{dt} = k_{hyd} \cdot (1 - H) \cdot \frac{\exp\left(\frac{T - 100}{15}\right)}{\rho}$$
*   **$k_{hyd}$**: Base hydration constant ($\approx 0.5 \text{ min}^{-1}$).
*   **Round Persistence**: Subsequent rounds ($n > 1$) start with the hydration level inherited from the previous session.

### 2.2. Thermodynamic Decay (Cooling)
$$\frac{dT}{dt} = -k_{cool} \cdot T + Q_{heat}$$
$$k_{cool} = 0.02 \cdot \frac{1}{C_p} \cdot \text{insulation factor}$$

*   **$C_p$**: Specific heat of vessel material (e.g., Zisha, Porcelain).
*   **Insulation**: Butter layer (Tibetan style) reduces $k_{cool}$ by 80%.

### 2.3. Component Extraction Kinetics
For each component $i$ (Catechin, Amino Acid, Caffeine, Pectin, Polysaccharide):
$$\frac{de_i}{dt} = k_{i,eff} \cdot (S_i - (e_{i,accum} + e_{i,curr}))$$
The effective rate constant $k_{i,eff}$ incorporates temperature, the hydration-driven geometry factor, water quality, and a session-based lag penalty:
$$k_{i,eff} = k_i \cdot \exp\left(sens_i \cdot \frac{T - 100}{30}\right) \cdot \sqrt{H} \cdot \omega_{water} \cdot \Gamma_{turbulence} \cdot lag(t_{\text{session}})$$
*   **$S_i$**: Maximum soluble mass in mg/g.
*   **$e_{i,accum}$**: Mass extracted in all previous infusions.
*   **$\omega_{water}$**: Water quality factor (see 2.6).
*   **$\Gamma_{turbulence}$**: Turbulence factor for active boiling (see 2.7).
*   **$lag(t_{session})$**: A 30-second ramp-up factor modeling fresh water contact.

#### Aromatic Capture and Amino Axes
*   **Aroma Extraction vs. Volatility**: Aromatic solubles are tracked with normal extraction kinetics, but a vessel-dependent volatility drain continuously removes a portion of the captured aromatics. Two axes are rendered in the TUI:
    *   `aroma_extraction_axis` (mg/min entering the liquor)
    *   `aroma_volatility_axis` (mg/min simultaneously lost to the air). Thin-lidded glass/gaiwan vessels and high turbulence drive this upward.
    *   The simulator now forecasts the next volatility step and, once that projection shows vapor losses climbing past roughly 40 % of the simultaneous extraction (≈60/40 split), it ends the pour—provided the captured aromatics are still above the human perception floor.
*   **Amino Body & Flow**: Amino acids expose a normalized body (`amino_depth_axis = e_{ami} / S_{ami}`) and a live slope (`amino_vibrancy_axis = de_{ami}/dt`). Together with the aroma axes they gate the priority-stop window so the cut happens only when the cup is still supple.
*   **Clarity Bias**: In gongfu/standard mode the simulator scales down pectin/polysaccharide kinetics using `clarity_bias = 0.55 + 0.45(1 - \text{saturation})`, suppressing heavy colloids so the liquor stays light even when multiple infusions are planned.

### 2.6. Water Quality Correction ($\omega_{water}$)
Solute solubility decreases as the solvent (water) is pre-saturated with minerals:
$$\omega_{water} = \max\left(0.7, 1.0 - \frac{TDS - 100}{2000}\right) \text{ if } TDS > 100$$

### 2.7. Active Boiling (Decoction) Model
When the vessel is actively heated ($T \ge T_{boil}$):
1.  **Turbulence Factor ($\Gamma_{turbulence}$)**: Forced convection increases mass transfer. $\Gamma = 3.0$ during active boil.
2.  **Vaporization (Concentration)**: Water volume ($V$) decreases, concentrating solutes.
    $$\frac{dV}{dt} = -0.5 \cdot \text{HeatLevel} \text{ (ml/min)}$$
3.  **Cell-Wall Rupture**: Structural integrity ($I$) degrades, releasing Lignin ($k_{lig}$).
    $$\frac{dI}{dt} = -0.15 \text{ (min}^{-1}\text{)}$$
    $$k_{lig} = 0.2 \cdot \left(1.0 - \frac{I}{0.3}\right) \cdot \Gamma \text{ (if } I < 0.3)$$

### 2.8. Termination Logic
The simulator employs two different stop triggers based on the brewing style:

#### A. Gongfu/Standard (Acceleration Surge)
Predicts the "balance point" where bitterness acceleration ($B''$) exceeds sweetness acceleration ($S''$):
$$\text{Stop if: } \frac{d^2B}{dt^2} > \frac{d^2S}{dt^2}$$
*   $B = 0.7 r_{cat} + 0.3 r_{caf}$
*   $S = 0.4 r_{ami} + 0.4 r_{pol} + 0.2 r_{pec}$
*   $r_i = e_{i,total} / S_i$ (Relative extraction of component $i$)

#### B. Decoction (Viscosity Index)
Stops when heavy solutes reach a target concentration representing a "thick" mouthfeel:
$$\text{Stop if: } \frac{(1.5e_{pec} + 1.0e_{pol} + 2.0e_{lig}) \cdot m_{leaf}}{V} \ge 3.0 \text{ mg/ml}$$

#### C. Aroma-Priority Stop Window
For gongfu/standard styles the simulator now looks at the aromatic and amino axes before the acceleration check:
*   **Capture Threshold**: Requires at least 62 % of the aromatic plateau with extraction still outpacing volatility by >0.003 mg/min.
*   **Preemptive Volatility Cut**: A forward-pass estimates the next volatility pulse and halts the pour whenever that projection shows volatilized aromatics rising to ~40 % of the simultaneous extraction (≈6:4 split), as long as the captured load is still above the sensory floor.
*   **Amino Guardrail**: Amino body must exceed 35 % with a vibrancy slope ≥0.015 mg/min so the liquor does not taste hollow.
*   **Clarity Ceiling**: Heavy colloids are monitored even outside decoction. If `(1.5e_{pec} + e_{pol} + 2e_{lig}) \cdot m_{leaf} / V` exceeds ~1.3 mg/ml, the simulator halts the infusion and allows extraction progress to dip accordingly, yielding the “clear, aromatic” profile that sits between British thickness and ultra-short Chinese rinses.
*   **Progress Override**: When either of the above style cues fire, extraction snaps to the override score so the dashboard reflects the intended airy-yet-fragrant mouthfeel.

### 2.9. Leaf Unfurling & Surface Area Geometry

The simulator treats leaf expansion as a mechanical response to hydration ($H$) and temperature ($T$):

1. **Unfurling Velocity ($u$):**
   $$\frac{du}{dt} = k_{unfurl} \cdot (H - u) \cdot \frac{T}{100}$$
   Where $u \in [0, 1]$ represents the "openness" of the leaf.

2. **Effective Surface Area Factor ($\eta$):**
   Unlike a simple linear model, extraction is constrained by the exposed surface area. We model this as:
   $$\eta = 0.2 + 0.8 \cdot u$$
   - When $u=0$ (tightly rolled), only 20% of the surface is in contact with water.
   - When $u=1$ (fully unfurled), 100% exposure is achieved.

3. **Coupled Extraction:**
   The effective rate constant $k_{i,eff}$ is directly scaled by $\eta$. This explains why rolled Oolongs or compressed Pu-erhs have a "slow start" regardless of temperature.

### 2.10. Vessel Geometry Factors ($\Xi$)

The physical shape and opening diameter of the vessel influence the extraction rate through agitation efficiency and surface-to-volume ratio:
$$k_{i,eff} \propto \Xi_{vessel}$$
- **$\Xi_{gaiwan} = 1.25$**: High opening-to-depth ratio allows for greater convective mixing and easier leaf expansion.
- **$\Xi_{teapot} = 1.00$**: Standard baseline for enclosed brewing.

### 2.11. Tea Distributor Algorithm

Long gongfu sessions can drag when one infusion stubbornly lags behind the rest, leading to uneven taste distribution across later pours. The **Tea Distributor Algorithm** continuously predicts how much soluble mass is still locked inside the leaves and automatically rebalances the extraction targets of the remaining infusions.

1. **Geometry-Aware Leaf Count**  
   The simulator reconstructs an approximate per-leaf mass using the configured leaf geometry (width & height), the tea profile's density, and a 0.3 mm thickness assumption:  
   $$m_{\text{leaf}} = \rho_{leaf} \times \frac{w \cdot h}{100} \times 0.03$$  
   The total leaf count $N$ is $m_{leaf\_mass} / m_{\text{leaf}}$. This feeds into a geometry factor  
   $$\Gamma_{geo} = 0.6 + 0.25 \frac{N}{N+4} + 0.15 \left(\tanh\frac{\max(w,h)}{45} + 0.5\right)$$

2. **Potential Concentration**  
   The remaining soluble load is computed from the component plateaus ($S_i$), tea mass, and total extraction level:  
   $$C_{pot} = \frac{\Gamma_{geo} \cdot ( \sum_i S_i ) \cdot m_{leaf} \cdot (1 - EY_{total})}{V_{inf}} \quad [\text{mg/ml}]$$  
   This value is exposed on the dashboard so the user can see how much flavor headroom is left.

3. **Drag Detection**  
   Expected kinetics are derived from the tea profile’s base velocity, hydration state, and the geometry factor. The real-time sweetness velocity ($dS/dt$) provides the measured rate. Their ratio is the **drag ratio**:
   $$\text{drag} = \text{clamp}\left(\frac{\dot{S}_{actual}}{\dot{S}_{expected}}, 0.4, 1.8\right)$$  
   Values below 1.0 indicate a lagging infusion.

4. **Dynamic Target Allocation**  
   For each infusion, the baseline target is the remaining extraction percentage divided by the number of infusions left. The distributor scales this slice by a composite factor that blends the potential concentration, tea-type velocity bias, leaf-length factor, and the drag ratio (clamped to \(0.4-1.6\times\)). The result becomes the new `target_extraction_for_cycle`, and the dashboard shows the per-infusion share hint.

This mechanism ensures that when one stage is slow or underperforming, the simulator automatically reallocates the remaining “flavor budget” so later infusions can still deliver meaningful cups.

---

## 3. Brewing Guidelines for Accurate Simulation

To achieve the best correlation between the physical simulation and your actual brew, follow these techniques:

### A. Controlled Water Flow
*   **Gentle Pouring**: For the most accurate prediction of hydration and leaf unfurling, **pour water slowly and steadily along the inner wall of the vessel**.
*   **Avoid Direct Impact**: Avoid pouring water directly onto the dry tea leaves at high velocity. This prevents uneven hydration (localized "shocking") and ensures the leaves unfurl according to the modeled laminar flow dynamics.
*   **Consistency**: A slow wall-pour maintains a more stable temperature profile during the initial contact phase, which the simulator assumes for its $lag(t_{session})$ calculations.

### B. Vessel Preheating
*   The thermodynamic model assumes the vessel has reached a partial thermal equilibrium. Briefly rinsing your teapot or gaiwan with hot water before adding leaves will align the real-world cooling rate with the simulator's $k_{cool}$ constant.

---

## 4. Implementation Details (libttak)

This project strictly adheres to the **libttak** memory safety model:
- **Lifetime Tracking**: Every allocation is associated with an explicit expiry.
- **Ownership**: The `tea_state_t` is managed as a `ttak_shared_t` resource with `TTAK_SHARED_LEVEL_3` security.
- **Memory Tree**: Real-time garbage collection is performed via `tt_autoclean_dirty_pointers`.
- **Safety Policy**: All data structure operations are isolated within a `ttak_owner_t` context.


### Build & Run
```bash
make
./bin/mobrew
```

### Dashboard Indicators
*   **Brew Progress**: Total EY relative to target.
*   **Leaf Hydration**: Real-time expansion status of the leaves.
*   **Astringency Δ**: Instantaneous rate of bitter compound release.
*   **Target Point**: The 100% mark where the balance of flavor is theoretically ideal.
*   **Aroma/Amino Axes**: Shows captured aroma mass plus extraction vs. volatility rates, and the amino body/flow pair that governs the aroma-priority stop window.
*   **Clarity Index**: Real-time heavy-colloid concentration (mg/ml) with a highlighted 1.3 threshold so you can see when the simulator is about to issue the “clear, aromatic” stop cue.

### Session Halt Guardrail
If the Tea Distributor determines the remaining flavor potential per infusion falls below ~0.35 % or the concentration potential drops under ~0.8 mg/ml, the dashboard surfaces a modal message, **“You can't steep more,”** and the session terminates—even if the user originally asked for seven or more infusions.
