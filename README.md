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
│       └── poppe.py       # Poppe method (ODE system)
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
from lib.cooling_tower.air import AirFlow
from lib.cooling_tower.water import WaterFlow
from lib.cooling_tower.poppe import PoppeSolver
from lib.cooling_tower.common import Q_, u

air = AirFlow(temp=Q_(25, u.degC), humidity=Q_(50, u.perc))
water = WaterFlow(temp=Q_(40, u.degC))
solver = PoppeSolver(air, water, water_temp_out=Q_(30, u.degC), lg_ratio=1.2)

results_df = solver.solve()
print(results_df.tail())
