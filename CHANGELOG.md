# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-05-07

### Added
- Unit-aware input normalization for core models via decorators in `units_descriptor.py`.
- Shared solver behavior in `SolverMixin`:
  - common solver initialization contract (`air_in`, `water_in`, `water_out`, `lg_ratio`, `C`, `n`)
  - shared target relation `target_me(lg_ratio)`
  - shared temperature search routine `estimate_temperatures(...)`.
- Extended test coverage for:
  - `AirFlow` psychrometric helpers and utility behavior
  - `WaterFlow` property trends and unit wrappers
  - `MerkelSolver` operating-point and temperature-estimation workflows
  - `PoppeSolver` profile, evaporation, and SNIP-inspired evaporation-order checks.
- Expanded `README.md` in English with:
  - architecture overview
  - method explanations (Merkel/Poppe)
  - numerical stability tricks
  - up-to-date usage examples and test commands.

### Changed
- `PoppeSolver.solve()` now returns a scalar Merkel-equivalent value (`float`).
- Poppe detailed outputs are exposed as solver instance state:
  - `profiles`
  - `evaporation`
  - `fog_carryover`
  - `total_loss`
  - `fog_force`.
- Relative humidity in Poppe profile output is stored in percent-like form (`air_rh_perc` = RH * 100).

### Removed
- Dedicated `find_temperatures_by_merkel(...)` implementation from `PoppeSolver`.
- Poppe temperature search is now unified through `estimate_temperatures(...)` from `SolverMixin`.

### Notes
- Core numerical path is float-first for performance/stability; quantity wrappers remain available through `calc_*` helper accessors.
- The `0.0012 * ΔT` evaporation relation is treated as an engineering reference check (SNIP guide context), not a strict equality constraint for detailed Poppe dynamics.
