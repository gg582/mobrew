# mobrew: Accurate, Simulated Tea Timer

[Demo](./demo.gif)

`mobrew` is not just a timer.
It’s a real-time kinetic simulator for the perfect cup of tea.

Most tea drinkers rely on static timers (e.g., "3 minutes for Green Tea"). However, tea extraction is a complex, non-linear physical process. `mobrew` uses a series of coupled differential equations to model how heat, leaf geometry, and water chemistry interact to produce the specific flavor profile in your cup.

## Why use a simulator instead of a timer?

Flavor is a race between different chemical compounds. Amino acids (sweetness/umami) often extract early and at lower temperatures, while catechins (bitterness/astringency) accelerate as the leaves unfurl and the temperature stays high.

`mobrew` tracks these variables in real-time, predicting the **"Flavor Equilibrium"**-the exact second where you maximize sweetness before bitterness takes over.

## Key Features

- **Leaf Unfurling Engine**: Models the mechanical opening of tea leaves to calculate the actual effective surface area.
- **Component-Specific Kinetics**: Tracks Catechins, Amino Acids, Caffeine, and Polysaccharides independently using Arrhenius-based scaling.
- **Thermodynamic Vessel Profiles**: Pre-configured cooling curves for Zisha (Purple Clay), Gaiwan, Porcelain, Glass, and British Teapots.
- **Water Chemistry Correction**: Adjusts extraction velocity based on TDS (Total Dissolved Solids) and mineral content.
- **Multi-Infusion Memory**: Tracks "leaf exhaustion" across multiple rounds (Gongfu style), so the 5th steep is as scientifically grounded as the 1st.

## The Secret to a Perfect Prediction

To get the most out of `mobrew`, your physical technique should match the model's assumptions:

1. **The "Wall Pour"**: Pour water slowly and steadily along the inner wall of your vessel. This ensures the laminar flow assumed by our hydration model and prevents "thermal shock" to the leaves.
2. **Preheat**: Always rinse your vessel with hot water first. This aligns the real-world thermal state with the simulator's starting constants.

## Installation

```bash
# Requires Ncurses and a C11 compiler
make
cp ./bin/mobrew /usr/local/bin
```

## Advanced Science

For a deep dive into the underlying math (Hydration Models, Mass Transfer Coefficients, and the `libttak` memory architecture), see [TECH.md](./TECH.md).
