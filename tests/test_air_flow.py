import numpy as np
import pint
from ..lib import *
import pytest as pt

def test_humidity_units():
    humidity = Q_(40.0, u.perc)
    assert humidity.to(u.ratio).magnitude == 0.4


def test_calc_wet_bulb_temperature():
    temp = Q_(30, u.degC)
    pressure = Q_(1, u.atm)
    humidity = Q_(40.0, u.perc)
    air = AirFlow(temperature=temp, humidity=humidity, pressure=pressure)
    t_wet = air.wet_bulb_temperature()
    assert t_wet == pt.approx(19.7, abs=0.5)

def test_omega_calculation():
    air = AirFlow(temperature=Q_(20, u.degC), humidity=Q_(50, u.perc), pressure=Q_(1, u.atm))
    w = air.omega()
    assert w == pt.approx(0.0073, abs=0.0001)

def test_enthalpy_increase_with_humidity():
    air_dry = AirFlow(temp=25.0, humidity=Q_(0, u.perc))
    air_wet = AirFlow(temp=25.0, rh=0.5)
    h1 = air_dry.wet_air_enthalpy()
    h2 = air_wet.wet_air_enthalpy()
    assert h2 > h1

def test_omega_at_zero_humidity():
    air = AirFlow(humidity=Q_(0, u.perc))
    assert air.omega() == pt.approx(0.0, abs=0.001)

@pt.mark.parametrize("temp_c, expected_kj", [
    (0, 0),       # Базовая точка
    (20, 20.12),  # 1.006 * 20
])
def test_dry_enthalpy_values(temp_c, expected_kj):
    air = AirFlow(temperature=Q_(temp_c, u.degC))
    h = air.dry_air_enthalpy()
    assert h/1000.0 == pt.approx(expected_kj, abs=0.01)

def test_wet_air_enthalpy():
    air = AirFlow(temperature=Q_(30, u.degC), humidity=Q_(60, u.perc))
    h_kj = air.wet_air_enthalpy()
    assert h_kj/1000.0 == pt.approx(71.0, abs=1.0)


def test_saturated_enthalpy_at_water_temp():
    air = AirFlow(pressure=Q_(101325, u.Pa))
    with air.temporary_set(temperature=Q_(35, u.degC)): 
        h_sat = air.saturated_air_enthalpy()/1000.0
    assert h_sat == pt.approx(128.8, abs=0.5)
    assert air.temp == 25.0

def test_saturated_omega_limit():
    air = AirFlow()
    with air.temporary_set(temperature=Q_(20, u.degC )):
        w_20 = air.saturated_omega()
    w_40 = air.saturated_omega(40.0)
    assert w_40 > w_20
    assert w_40 == pt.approx(0.049, abs=0.002)

def test_merkel_driving_force():
    t_air = Q_(25, u.degC)
    rh_air = Q_(50, u.perc)
    t_water = Q_(35, u.degC) 
    air = AirFlow(temperature=t_air, humidity=rh_air)
    h_air = air.wet_air_enthalpy()/1000.0
    h_sat_water = air.saturated_air_enthalpy(t_water)/1000.0
    assert h_sat_water > h_air
    driving_force = h_sat_water - h_air
    assert driving_force > 50

def test_poppe_omega_difference():
    air = AirFlow(temperature=Q_(25, u.degC), humidity=Q_(50, u.perc))
    t_water = Q_(35, u.degC)  
    w_air = air.omega()
    w_sat = air.saturated_omega(t_water)
    assert w_sat > w_air
    delta_w = w_sat - w_air
    assert delta_w > 0

def test_lewis_factor_range():
    """Число Льюиса должно быть в разумных пределах для атмосферного воздуха."""
    air = AirFlow(temp=Q_(25, u.degC), humidity=Q_(60, u.perc))
    le = air.lewis_factor()
    
    # Физически обоснованный диапазон для градирен
    assert 0.8 <= le <= 1.0

def test_calc_wet_bulb_returns_quantity():
    air = AirFlow(temp=30.0, rh=0.4, press=101325.0)
    q = air.calc_wet_bulb_temperature()
    assert q.units == u.degC
    plain = air.wet_bulb_temperature()
    assert q.magnitude == pt.approx(float(getattr(plain, "magnitude", plain)), abs=1e-6)

def test_vapor_enthalpy_linear_in_temp():
    air = AirFlow()
    t_c = 30.0
    hv = air.vapor_enthalpy(t_c)
    assert hv == pt.approx(air.EVW + air.CPV * t_c, rel=1e-6)

def test_calc_temperature_from_h_w_round_trip():
    air = AirFlow(temp=25.0, rh=0.5, press=101325.0)
    w = air.omega()
    h = psychrolib.GetMoistAirEnthalpy(25.0, w)
    t_back = air.temperature_from_h_w(h_j_kg=h, omega_kg_kg=w)
    assert t_back == pt.approx(25.0, abs=0.05)

def test_calc_temperature_from_h_w_with_units_flag():
    air = AirFlow(temp=20.0, rh=0.5, press=101325.0)
    w = air.omega()
    h = psychrolib.GetMoistAirEnthalpy(20.0, w)
    q = air.calc_temperature_from_h_w(h_j_kg=h, omega_kg_kg=w)
    assert hasattr(q, "magnitude")
    assert q.magnitude == pt.approx(20.0, abs=0.05)

def test_temporary_set_restores_and_rejects_unknown():
    air = AirFlow(temp=20.0, rh=0.5)
    orig_t = air.temp
    with air.temporary_set(rh=0.9):
        assert air.rh == pt.approx(0.9)
    assert air.rh == pt.approx(0.5)
    assert air.temp == orig_t
    with pt.raises(AttributeError):
        with air.temporary_set(not_a_field=1):
            pass