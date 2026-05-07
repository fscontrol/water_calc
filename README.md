# Cooling Tower Calculation Library

A Python library for cooling tower thermal calculations using two approaches:
- **Merkel method** (integral model)
- **Poppe method** (differential model with humidity profile and fog handling)

The codebase is designed around a practical rule:
- **Core state and numerics run on plain SI floats** (stable and fast)
- **Unit decorators allow convenient engineering inputs** (`pint` quantities)

---

## Project Layout

```text
.
├── lib/
│   ├── __init__.py
│   └── cooling_tower/
│       ├── __init__.py           # Public exports
│       ├── air.py                # Moist air / psychrometric helpers
│       ├── water.py              # Water/brine properties
│       ├── merkel.py             # Merkel solver
│       ├── poppe.py              # Poppe solver
│       ├── mixins.py             # Shared solver and utility mixins
│       ├── units_descriptor.py   # Unit descriptors + decorators
│       └── common.py             # Unit registry and common imports
├── tests/                        # Pytest suite
├── merkel.ipynb                  # Merkel examples
├── poppe.ipynb                   # Poppe examples
├── pixi.toml                     # Environment definition
└── README.md
```

---

## Main Concepts

### 1) Float-first internal model

All heavy computations (ODE integration, root finding, correlations) are expected to run with plain floats in SI-compatible units:
- temperature in `degC`
- pressure in `Pa`
- humidity ratio in `kg/kg`
- enthalpy in `J/kg`

This avoids performance and stability issues from carrying `Quantity` objects through every numeric call.

### 2) Unit-aware API layer

The project uses two decorators to keep APIs convenient while preserving float internals:
- `@cleans(...)` and `@cleans_simple(...)`: convert incoming `pint` values to magnitudes in target units
- `@returns(unit)`: declares return unit metadata for `calc_*` wrappers

`UnitMagicMixin` provides lazy wrappers such as:
- `obj.calc_density()`
- `obj.calc_wet_bulb_temperature()`

These wrappers return `pint.Quantity` without forcing the main computation path to use quantities.

### 3) Shared solver behavior in `SolverMixin`

`SolverMixin` centralizes shared solver state and helper logic:
- common constructor signature: `air_in, water_in, water_out, lg_ratio, C, n`
- target Merkel relation: `target_me(lg_ratio) = C * lg_ratio**(-n)`
- generic temperature estimation routine: `estimate_temperatures(...)`

This allows both `MerkelSolver` and `PoppeSolver` to reuse the same temperature-search workflow.

---

## Physics Models

## Merkel Solver (`merkel.py`)

- Integrates classical Merkel expression across water temperature range.
- Uses a driving force based on saturated air enthalpy minus current air enthalpy line.
- Supports:
  - direct Merkel number computation via `solve(...)`
  - operating `L/G` search via `find_operating_lg()`
  - inlet/outlet temperature estimation via shared `estimate_temperatures(...)`

## Poppe Solver (`poppe.py`)

- Solves coupled ODEs for air enthalpy and humidity ratio along water-temperature coordinate.
- Handles two regimes:
  - unsaturated air
  - fog/saturated branch
- Produces both scalar and profile outputs:
  - `solve()` returns a scalar Merkel-equivalent number (`float`)
  - side effects on solver instance:
    - `profiles` (`DataFrame` with fields like `air_temp_c`, `air_omega_kg_kg`, `zone`, etc.)
    - `evaporation`, `fog_carryover`, `total_loss`, `fog`

---

## Practical Numerical Tricks Used

The implementation includes several pragmatic engineering choices to keep solvers robust:

1. **Flooring small denominators**
   - Expressions like `max(denom, 1e-7)` prevent division blowups near equilibrium.

2. **Bounded integration step**
   - `solve_ivp(..., max_step=0.5)` improves stability for stiff or rapidly changing regions.

3. **Stiff solver for Poppe ODEs**
   - `method='Radau'` is selected for robustness in nonlinear coupled equations.

4. **Bracketed root finding**
   - `brentq` is used where sign change is guaranteed/checked.
   - This is safer than open methods when residuals are expensive and sometimes irregular.

5. **Temporary state mutation with rollback**
   - `temporary_set(...)` context manager allows trying candidate temperatures safely and restoring state automatically.

6. **Two-layer API design**
   - fast float path for internals
   - `calc_*` wrappers for unit-aware external use

---

## Installation and Environment (Pixi)

Using [Pixi](https://pixi.sh):

```bash
# clone project, then
cd /path/to/water_calc
pixi install

# run tests
pixi run pytest

# open notebooks
pixi run jupyter lab
```

---

## Quick Usage Examples

### Merkel

```python
from lib import AirFlow, WaterFlow, MerkelSolver, Q_, u

air = AirFlow(temperature=Q_(25, u.degC), humidity=Q_(50, u.perc))
water_in = WaterFlow(temperature=Q_(40, u.degC))
water_out = WaterFlow(temperature=Q_(30, u.degC))

solver = MerkelSolver(air_in=air, water_in=water_in, water_out=water_out, lg_ratio=1.0, C=1.0, n=0.6)
me = solver.solve()

# shared temperature search
t_in, t_out, err = solver.estimate_temperatures(lg_ratio=1.0, delta_t=Q_(10, u.delta_degC))
```

### Poppe

```python
from lib import AirFlow, WaterFlow, PoppeSolver, Q_, u

air = AirFlow(temperature=Q_(25, u.degC), humidity=Q_(60, u.perc))
water_in = WaterFlow(temperature=Q_(40, u.degC))
water_out = WaterFlow(temperature=Q_(30, u.degC))

solver = PoppeSolver(air_in=air, water_in=water_in, water_out=water_out, lg_ratio=1.2)
me = solver.solve()  # float

# profile and losses are stored on the solver object
df = solver.profiles
evap = solver.evaporation
fog = solver.fog_carryover
total = solver.total_loss
```

---

## Testing

Run all tests:

```bash
pixi run pytest
```

Run specific suites:

```bash
pixi run pytest tests/test_merkel.py -v
pixi run pytest tests/test_poppe.py -v
pixi run pytest tests/test_air_flow.py tests/test_water_flow.py -v
```

---

## Notes and Conventions

- `humidity` is treated as relative humidity ratio (`0..1`), with convenience unit `u.perc`.
- `WaterFlow.salinity` is represented as `g/kg` equivalent scale in current correlations.
- Poppe profile RH column is kept in percent-like naming (`air_rh_perc`) for reporting.
- When extending the library, prefer:
  1. keep internal solvers float-based
  2. normalize inputs at API boundaries (`@cleans*`)
  3. expose unit-aware views through `calc_*` wrappers

This pattern keeps the code numerically stable while still ergonomic in notebooks.
