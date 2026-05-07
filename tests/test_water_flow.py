import pint
from ..lib import *
import pytest as pt

def test_brine_properties():
    sea_water = WaterFlow(temperature=Q_(30, u.degC), tds=Q_(35, u.g / u.kg))
    assert sea_water.density() == pt.approx(1022, abs=5)
    fresh_water = WaterFlow(temperature=Q_(30, u.degC), tds=Q_(0, u.g / u.kg))
    assert sea_water.viscosity() > fresh_water.viscosity()
    assert sea_water.latent_heat() < fresh_water.latent_heat()

def test_latent_heat_at_high_salinity():
    brine = WaterFlow(temperature=Q_(40, u.degC), tds=Q_(100, u.g / u.kg))
    lv = brine.calc_latent_heat().to("kJ/kg")
    assert lv.magnitude < 2400

def test_specific_heat_decreases_with_salinity():
    t = 35.0
    fresh = WaterFlow(temp=t, salinity=0.0)
    brine = WaterFlow(temp=t, salinity=50.0)
    assert brine.specific_heat() < fresh.specific_heat()

def test_calc_density_wraps_quantity():
    w = WaterFlow(temp=25.0, salinity=0.0)
    rho = w.calc_density()
    assert rho.units == u.kg / u.m**3
    assert rho.magnitude == pt.approx(997.0, abs=2)

def test_viscosity_positive():
    w = WaterFlow(temp=20.0, salinity=35.0)
    assert w.viscosity() > 0
    assert w.calc_viscosity().magnitude > 0