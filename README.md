# Cooling Tower Calculation Library

Python library for thermal design checks of wet cooling towers using:

- **Merkel** — classical integral (enthalpy driving force) model  
- **Poppe** — differential model along water temperature with humidity evolution and a fog / supersaturation branch  

Design principle:

- **Internals:** plain SI floats (fast, stable with `scipy` / `psychrolib`)  
- **API:** optional `pint` inputs via descriptors and `@cleans` / `@cleans_simple` / `@returns`; `calc_*` helpers on mixins return `Quantity` objects  

---

## Layout

```text
.
├── lib/
│   └── cooling_tower/
│       ├── air.py              # Moist air, psychrometrics, moist-air density
│       ├── water.py            # Water / brine properties, optional flow / evaporation hooks
│       ├── merkel.py           # Merkel solver
│       ├── poppe.py            # Poppe ODE solver (solve_ivp)
│       ├── mixins.py           # SolverMixin, TemporarySetMixin, UnitMagicMixin
│       ├── units_descriptor.py # Unit descriptors + decorators
│       └── common.py           # pint registry, shared imports
├── tests/
├── merkel.ipynb, poppe.ipynb
├── pixi.toml
├── CHANGELOG.md
└── README.md
```

Public imports: `from lib import AirFlow, WaterFlow, MerkelSolver, PoppeSolver, Q_, u` (see `lib/cooling_tower/__init__.py`).

---

## Core ideas

### Float-first numerics

Use **°C**, **Pa**, humidity ratio **kg/kg dry air**, enthalpy **J/kg** inside solvers.  
Avoid passing `pint.Quantity` into `psychrolib` or tight loops.

### Unit-aware boundary

- `@cleans` / `@cleans_simple` — normalize arguments to magnitudes in declared units  
- `@returns(unit)` — metadata for `UnitMagicMixin` so `obj.calc_*()` wraps plain-float methods as `Quantity`  

### Shared solver API (`SolverMixin`)

Constructor (positional or keywords):

`AirFlow air_in`, `WaterFlow water_in`, `WaterFlow water_out`, `lg_ratio`, `C`, `n`

- **`lg_ratio`** — liquid-to-gas mass ratio \(\dot m_w / \dot m_{da}\) (water mass per **dry** air mass), dimensionless in the code.  
- **`target_me(lg)`** — manufacturer-style characteristic: `C * lg**(-n)` (Merkel number target vs. operating \(L/G\)).  
- **`find_operating_lg()`** — finds \(L/G\) such that `solve(lg) ≈ target_me(lg)` (bracket + `brentq`).  
- **`estimate_temperatures(lg_ratio=..., delta_t=..., target_me=..., return_units=False)`** — finds cold / hot water temperatures with fixed approach \(\Delta T = T_{in} - T_{out}\) so that `solve` matches the target Merkel.  
  - Lower search limit starts at `Twb + 2 °C` and is shifted upward in small steps while `solve` raises `ValueError` (e.g. stiff Poppe / psychrometric domain), until a valid `f(a)` is obtained.  
  - Returns `(t_in, t_out)` as floats, or `(Q_, Q_)` if `return_units=True`.  

`TemporarySetMixin.temporary_set(**kwargs)` temporarily overrides attributes (e.g. water temperatures) and restores them after the block — used inside `estimate_temperatures`.

---

## Merkel (`merkel.py`)

- Integrates \(\int_{T_{out}}^{T_{in}} \frac{c_{p,w}(T_w)}{h_{sat}(T_w) - h_a(T_w)}\,\mathrm d T_w\) with a floor on the denominator (`max(..., 0.1)`).  
- **`solve(lg_ratio=None)`** — Merkel number for current states.  
- **`air_enthalpy_at_tw(t_w, lg_ratio=None)`** — air enthalpy along the operating line.  
- Fails if \(T_{out}\) is below wet-bulb (non-physical outlet).

---

## Poppe (`poppe.py`)

- State along water temperature \(T_w\) (counter-flow: integrate from **`water_out.temp` → `water_in.temp`**).  
- **`poppe_system`** — RHS for \(d h_a / dT_w = (L/G)\, c_{p,w}\) and \(d\omega/dT_w\) from Poppe’s coupling; unsaturated branch vs. fog branch when \(\omega \ge \omega_{sat}(T_{air})\).  
- **`solve(lg_ratio=None, steps=100)`** — `solve_ivp(..., method='Radau', max_step=0.5, rtol/atol as in code)` on a uniform \(T_w\) grid of `steps` points.  
- **Return value:** Merkel-equivalent integral (trapezoid over the same profile) as `float`.  
- **Side effects on the instance:**  
  - **`profiles`** — `DataFrame`: `water_temp_c`, `air_temp_c`, `air_rh_perc`, `air_omega_kg_kg`, `air_enthalpy_j_kg`, `fog_kg_kg`, `zone`  
  - **`evaporation`**, **`fog_carryover`**, **`total_loss`** — per unit \(L/G\) convention used in code  
  - **`fog_force`** — fraction of profile points classified as fog (×100 for a percentage-like number)  
- **`decode_results`** — recovers dry-bulb and RH from \((h,\omega)\) at fixed pressure.

Numerical safeguards: `max(denom, 1e-7)` in the unsaturated \(\mathrm d\omega/\mathrm dT_w\) expression; stiff ODE integrator; bracketed root finding where applicable.

---

## Air (`air.py`)

Psychrometrics via **`psychrolib`** (SI): wet bulb, \(\omega\), enthalpies, saturation, Lewis factor, \(\mathrm d\omega_{sat}/\mathrm dT\), etc.

- **`temperature_from_h_w`** — approximate \(T(h,\omega)\) from ideal gas / latent heat line (used inside Poppe for dry-bulb state).  
- **`density(temp, rh, press)`** — **`GetMoistAirDensity(Tdb, ω, P)`** with \(\omega\) from `GetHumRatioFromRelHum` (second argument to psychrolib must be humidity ratio, not RH).  
- **`volume_to_dry_mass(volume, ...)`** — \(\dot m_{da} \approx \rho_{moist}\,\dot V / (1+\omega)\) for consistent \(L/G\) from volumetric fan data.  
- Optional **`flow`** on `AirFlow` for bookkeeping (not required by solvers).

**Note:** `@returns(u.g/u.kg)` on methods that actually return **kg/kg** from psychrolib is only relevant for `calc_*` display; internal values remain kg/kg floats.

---

## Water (`water.py`)

Brine / water: **`density`**, **`viscosity`**, **`latent_heat`**, **`specific_heat`** vs. `temp` and `salinity` (g/kg scale in correlations).  
Optional **`flow`** and **`apply_evaporation()`** for simple mass / salinity updates (evaporation ratio must be set on the instance by your workflow).

---

## Mass \(L/G\) from a fan (volumetric flow)

Fans are often rated in m³/s; tower \(L/G\) is mass-based. At air inlet \((T_{db}, \varphi, P)\):

1. \(\omega =\) `GetHumRatioFromRelHum(Tdb, φ, P)`  
2. \(\rho_{moist} =\) `GetMoistAirDensity(Tdb, ω, P)`  
3. \(\dot m_{moist} = \rho_{moist}\,\dot V\)  
4. \(\dot m_{da} = \dot m_{moist} / (1+\omega)\)  
5. \(L/G = \dot m_w / \dot m_{da}\)  

`AirFlow.volume_to_dry_mass` implements step 4 given \(\dot V\) and state.

---

## Environment (Pixi)

```bash
cd /path/to/water_calc
pixi install
pixi run pytest
pixi run jupyter lab   # notebooks
```

Dependencies: `numpy`, `scipy`, `pandas`, `pint`, `psychrolib`, `pytest`, `jupyter`, `matplotlib`; see `pixi.toml`.

---

## Examples

### Merkel

```python
from lib import AirFlow, WaterFlow, MerkelSolver, Q_, u

air = AirFlow(temp=Q_(25, u.degC), rh=0.5, press=101325.0)
water_in = WaterFlow(temp=40.0)
water_out = WaterFlow(temp=30.0)

solver = MerkelSolver(air, water_in, water_out, lg_ratio=1.0, C=1.0, n=0.6)
me = solver.solve()

t_in, t_out = solver.estimate_temperatures(lg_ratio=1.0, delta_t=Q_(10, u.delta_degC))
```

### Poppe

```python
from lib import AirFlow, WaterFlow, PoppeSolver, Q_, u

air = AirFlow(temp=Q_(25, u.degC), humidity=Q_(60, u.perc))
water_in = WaterFlow(temp=40.0)
water_out = WaterFlow(temp=30.0)

solver = PoppeSolver(air, water_in, water_out, lg_ratio=1.2, C=1.0, n=0.6)
me = solver.solve(steps=100)

prof = solver.profiles          # pandas DataFrame
evap = solver.evaporation
```

---

## Tests

```bash
pixi run pytest
pixi run pytest tests/test_poppe.py -v
pixi run pytest tests/test_merkel.py -v
```

---

## Conventions

- Relative humidity: **0–1** float internally; `humidity=Q_(60, u.perc)` is supported on `AirFlow`.  
- **`WaterFlow`**: constructor uses `temp` / `salinity` (floats after `@cleans_simple`); descriptors `temperature` / `tds` exist for quantity-style assignment in callers that use them.  
- Extending: keep solver cores float-based; convert at boundaries; add `calc_*` via `@returns` when exposing engineering units.
