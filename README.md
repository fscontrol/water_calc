# Engineering Suite: Cooling Tower Thermodynamics & Water Chemistry
##  Project Overview

This project is a comprehensive tool for modeling industrial cooling towers, circulating water balances, and water treatment systems.

### Current Progress
Currently implemented:
* **AirFlow:** Moist air property calculations (enthalpy, humidity ratio, wet-bulb temperature) using `psychrolib`.
* **WaterFlow:** Physical properties of water and brines (density, viscosity, heat capacity) based on salinity.
* **Merkel Method:** Classical Merkel Number ($Me$) calculations and cooling performance forecasting.

### Development Roadmap
1.  **Poppe Method:** Transition from the simplified Merkel approach to a rigorous model (accounting for water mass loss and Lewis factor $Le \neq 1$).
2.  **Chemistry with REAKTORO:** Modeling chemical equilibria in circulating water. Predicting saturation indices based on makeup water composition.
3.  **Reverse Osmosis (RO):** A module for membrane performance normalization and makeup water desalination calculations.

---

## Tech Stack
* **Python 3.10+**
* **Pint:** Strict unit handling (SI and IP).
* **SciPy/NumPy:** Numerical integration and solvers.
* **Reaktoro:** Thermodynamic modeling of aqueous systems.
* **Psychrolib:** Standardized psychrometric calculations.
* **PIXI** package and virtual environment manager for Python
