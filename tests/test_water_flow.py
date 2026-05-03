import pint
from ..lib import *
import pytest as pt

def test_brine_properties():
    sea_water = WaterFlow(temp=Q_(30, u.degC), salinity=Q_(35, u.g/u.kg))
    assert sea_water.density().magnitude == pt.approx(1022, abs=5)
    fresh_water = WaterFlow(temp=Q_(30, u.degC), salinity=Q_(0, u.g/u.kg))
    assert sea_water.viscosity() > fresh_water.viscosity()
    assert sea_water.latent_heat() < fresh_water.latent_heat()

def test_latent_heat_at_high_salinity():
    brine = WaterFlow(temp=Q_(40, u.degC), salinity=Q_(100, u.g/u.kg))
    lv = brine.latent_heat().to('kJ/kg')
    assert lv.magnitude < 2400