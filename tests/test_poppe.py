import pytest
import numpy as np
from ..lib import *

# Пособие к СНиП 2.04.02-84 (градирни), п. 3.17, формула (22):
#   G_и = k_и · Δt · G_ж  →  удельные потери G_и/G_ж = k_и · Δt .
# k_и — по графику (черт. 61) от t воздуха по сухому термометру; для t_сух ≈ 20 °С
# в практике часто используют ориентир k_и ≈ 1,2·10⁻³ 1/°C (не норма точного числа).
def snip_evaporation_ratio(delta_t_celsius: float, k_i_per_deg_c: float = 0.0012) -> float:
    """Доля испарённой воды на долю циркуляции (кг/кг ≈ м³/м³), по упрощённому ориентиру."""
    return k_i_per_deg_c * float(delta_t_celsius)


def test_poppe_solver_execution():
    """Check if the solver runs and returns a non-empty DataFrame."""
    air_in = AirFlow(temperature=Q_(25, u.degC), humidity=Q_(60, u.perc))
    water_in = WaterFlow(temperature=Q_(40, u.degC))
    water_out = WaterFlow(temperature=Q_(30, u.degC))
    lg = 1.0

    solver = PoppeSolver(air_in, water_in, water_out, lg)
    df = solver.solve()

    assert not df.empty
    assert 'air_enthalpy_j_kg' in df.columns
    assert 'water_temp_c' in df.columns

def test_poppe_physics_trends():
    """Air enthalpy should increase as it moves through the tower."""
    air_in = AirFlow(temperature=Q_(20, u.degC), humidity=Q_(50, u.perc))
    water_in = WaterFlow(temperature=Q_(45, u.degC))
    water_out = WaterFlow(temperature=Q_(32, u.degC))

    solver = PoppeSolver(air_in, water_in, water_out, lg_ratio=1.0)
    df = solver.solve()

    h_start = df['air_enthalpy_j_kg'].iloc[0]
    h_end = df['air_enthalpy_j_kg'].iloc[-1]

    # In a cooling tower, air absorbs heat
    assert h_end > h_start

def test_poppe_mass_balance():
    """Evaporation should be positive in standard conditions."""
    air_in = AirFlow(temperature=Q_(30, u.degC), humidity=Q_(30, u.perc))
    water_in = WaterFlow(temperature=Q_(45, u.degC))
    water_out = WaterFlow(temperature=Q_(35, u.degC))
    solver = PoppeSolver(air_in, water_in, water_out, lg_ratio=1.1)
    df = solver.solve()

    # Check attributes stored in metadata
    assert df.attrs['evaporation_kg_kg'] > 0
    assert df.attrs['total_water_loss_kg_kg'] >= df.attrs['evaporation_kg_kg']


def test_poppe_evaporation_vs_snip_order_of_magnitude_at_20c():
    """Сравнение удельного испарения Поппе с ориентиром k_и·Δt (пособие к СНиП, ф. 22)."""
    t_out = 30.0
    t_in = 40.0
    delta_t = t_in - t_out
    air_in = AirFlow(temperature=Q_(20, u.degC), humidity=Q_(60, u.perc))
    water_in = WaterFlow(temp=t_in)
    water_out = WaterFlow(temp=t_out)
    df = PoppeSolver(air_in, water_in, water_out, lg_ratio=1.0).solve()
    e_poppe = float(df.attrs["evaporation_kg_kg"])
    e_snip = snip_evaporation_ratio(delta_t)
    assert e_poppe > 0
    # Детальная модель не обязана совпадать с графиком k_и; ожидаем тот же порядок
    assert e_poppe == pytest.approx(e_snip, rel=0.35)


def test_poppe_evaporation_increases_with_water_delta_t():
    """При прочих равных больший перепад воды даёт не меньшее испарение (Поппе)."""
    air_in = AirFlow(temperature=Q_(20, u.degC), humidity=Q_(50, u.perc))
    t_out = 28.0
    df_small = PoppeSolver(
        air_in, WaterFlow(temp=38.0), WaterFlow(temp=t_out), 1.0
    ).solve()
    df_large = PoppeSolver(
        air_in, WaterFlow(temp=45.0), WaterFlow(temp=t_out), 1.0
    ).solve()
    assert df_large.attrs["evaporation_kg_kg"] >= df_small.attrs["evaporation_kg_kg"]


def test_snip_evaporation_ratio_numeric_examples():
    assert snip_evaporation_ratio(10.0) == pytest.approx(0.012)
    assert snip_evaporation_ratio(12.0, k_i_per_deg_c=0.0012) == pytest.approx(0.0144)
