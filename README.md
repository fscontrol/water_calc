# Cooling Tower Calculation Library

A Python library for thermodynamic calculations of cooling towers using **Merkel** and **Poppe** methods. 
The library supports complex calculations including fog formation, evaporation losses, and variable water salinity.

## Project Structure
```text
.
├── lib/
│   └── cooling_tower/
│       ├── __init__.py    # Exports all modules
│       ├── air.py         # Air properties and psychrometrics
│       ├── common.py      # Units (Pint) and Shared constants
│       ├── merkel.py      # Merkel method implementation
│       ├── water.py       # Water/Brine properties
│       ├── poppe.py       # Poppe method (ODE system)
|       └── __init__.py
|    └──__init__.py
├── tests/                 # Pytest suite
├── merkel.ipynb           # Examples for Merkel method
├── poppe.ipynb            # Examples for Poppe method
├── pixi.toml              # Environment management
└── README.md
```

## Features
- **Psychrometrics:** Powered by `psychrolib` (SI units).
- **Unit Safety:** All inputs use `pint` for dimensional analysis.
- **Poppe Method:** Solves a system of ODEs to account for water evaporation and fogging.
- **Brine Support:** Includes salinity effects on density, viscosity, and latent heat.

## Quick Start
Using [Pixi](https://pixi.sh):

```bash
# install pixi
#for LINUX  (developed and tested under LINUX)
curl -fsSL https://pixi.sh/install.sh | bash
#for Windows
iwr -useb https://pixi.sh/install.ps1 | iex
#install dependencies
cd /path/to/downloaded/project
pixi install
# Run tests
pixi run pytest
# Launch notebooks
pixi run jupyter lab
```

## Usage Example
```python
from lib import *

air = AirFlow(temp=Q_(25, u.degC), humidity=Q_(50, u.perc))
water = WaterFlow(temp=Q_(40, u.degC))
solver = PoppeSolver(air, water, water_temp_out=Q_(30, u.degC), lg_ratio=1.2)

results_df = solver.solve()
print(results_df.tail())
#df structure
#            "water_temp_c":      temperatures of water alongside filling
#            "air_temp_c":        temperatures of air
#            "air_rh_perc":       humidity of air 
#            "air_omega_kg_kg":   vapor/dry air relation
#            "air_enthalpy_j_kg": air enthalpy
#            "fog_kg_kg":         kg of fog per kg of dry air
#            "zone":              zones (fog or unsatureted, string)
# also df has attributes can be accessed as dict
# df.attrs["total_water_loss_kg_kg"] 
# --- "merkel_number"
# --- "fog_carryover_kg_kg"
# --- "evaporation_kg_kg"
results = PopperSolver.process_results(df)
#will return
#    dict(
#            evaporation = df.attrs["total_water_loss_kg_kg"],
#            me = df.attrs["merkel_number"],
#            fog = df.attrs["fog_carryover_kg_kg"],
#            fog_force = float((df["zone"] == "fog").sum()/df["zone"].size*100.0))
