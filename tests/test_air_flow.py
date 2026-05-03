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
    air = AirFlow(temp=temp, humidity=humidity, pressure=pressure)
    t_wet = air.wet_bulb_temperature()
    assert t_wet.units == u.degC
    assert t_wet.magnitude == pt.approx(19.7, abs=0.5)

def test_omega_calculation():
    air = AirFlow(temp=Q_(20, u.degC), humidity=Q_(50, u.perc), pressure=Q_(1, u.atm))
    w = air.omega()
    assert w.units == u.ratio
    assert w.magnitude == pt.approx(0.0073, abs=0.0001)

def test_enthalpy_increase_with_humidity():
    t = Q_(25, u.degC)
    air_dry = AirFlow(temp=t, humidity=Q_(0, u.perc))
    air_wet = AirFlow(temp=t, humidity=Q_(50, u.perc))
    h1 = air_dry.wet_air_enthalpy()
    h2 = air_wet.wet_air_enthalpy()
    assert h2 > h1

def test_omega_at_zero_humidity():
    air = AirFlow(humidity=Q_(0, u.perc))
    assert air.omega().magnitude == pt.approx(0.0, abs=0.001)

@pt.mark.parametrize("temp_c, expected_kj", [
    (0, 0),       # Базовая точка
    (20, 20.12),  # 1.006 * 20
])
def test_dry_enthalpy_values(temp_c, expected_kj):
    air = AirFlow(temp=Q_(temp_c, u.degC))
    h = air.dry_air_enthalpy().to(u.kJ / u.kg)
    assert h.magnitude == pt.approx(expected_kj, abs=0.01)

def test_units_conversion_in_tests():
    # Проверяем, что класс корректно возвращает кДж/кг при запросе
    air = AirFlow(temp=Q_(30, u.degC), humidity=Q_(60, u.perc))
    h_kj = air.wet_air_enthalpy().to('kJ/kg')
    assert h_kj.magnitude == pt.approx(71.0, abs=1.0)


def test_saturated_enthalpy_at_water_temp():
    air = AirFlow(pressure=Q_(101325, u.Pa))
    t_water = Q_(35, u.degC)
    h_sat = air.saturated_air_enthalpy(t_water).to('kJ/kg')
    assert h_sat.magnitude == pt.approx(128.8, abs=0.5)

def test_saturated_omega_limit():
    air = AirFlow()
    w_20 = air.saturated_omega(Q_(20, u.degC))
    w_40 = air.saturated_omega(Q_(40, u.degC))
    assert w_40 > w_20
    assert w_40.magnitude == pt.approx(0.049, abs=0.002)

def test_merkel_driving_force():
    t_air = Q_(25, u.degC)
    rh_air = Q_(50, u.perc)
    t_water = Q_(35, u.degC) 
    air = AirFlow(temp=t_air, humidity=rh_air)
    h_air = air.wet_air_enthalpy().to('kJ/kg')
    h_sat_water = air.saturated_air_enthalpy(t_water).to('kJ/kg')
    assert h_sat_water > h_air
    driving_force = h_sat_water - h_air
    assert driving_force.magnitude > 50

def test_poppe_omega_difference():
    air = AirFlow(temp=Q_(25, u.degC), humidity=Q_(50, u.perc))
    t_water = Q_(35, u.degC)  
    w_air = air.omega()
    w_sat = air.saturated_omega(t_water)
    assert w_sat > w_air
    delta_w = w_sat - w_air
    assert delta_w.magnitude > 0

def test_lewis_factor_range():
    """Число Льюиса должно быть в разумных пределах для атмосферного воздуха."""
    air = AirFlow(temp=Q_(25, u.degC), humidity=Q_(60, u.perc))
    le = air.lewis_factor()
    
    # Физически обоснованный диапазон для градирен
    assert 0.8 <= le <= 1.0