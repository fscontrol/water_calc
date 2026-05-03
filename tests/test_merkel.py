import pint
from ..lib import *
import pytest as pt

def test_merkel_logic():
    air = AirFlow(temp=Q_(25, u.degC), humidity=Q_(50, u.perc))
    water = WaterFlow(temp=Q_(40, u.degC)) # Вход 40
    t_out = Q_(30, u.degC)  
    solver = MerkelSolver(1.0, 0.6, air, water, t_out)
    me_val = solver.solve_me(lg_ratio=1.0)
    assert me_val > 0
    solver_colder = MerkelSolver(1.0, 0.6, air, water, Q_(28, u.degC))
    me_colder = solver_colder.solve_me(lg_ratio=1.0)
    
    assert me_colder > me_val

def test_merkel_standard_case():
    air_in = AirFlow(temp=Q_(25, u.degC), humidity=Q_(50, u.perc))
    water_in = WaterFlow(temp=Q_(40, u.degC))
    t_out = Q_(28, u.degC)
    solver = MerkelSolver(1.0, 0.6, air_in, water_in, t_out)
    me_val = solver.solve_me(lg_ratio=1.0)
    assert 1.0 < me_val < 1.5
    print(f"\nРассчитанное число Меркеля: {me_val:.4f}")

def test_merkel_limit_approach():
    t_dry = Q_(25, u.degC)
    rh = Q_(60, u.perc)
    air = AirFlow(temp=t_dry, humidity=rh)
    t_wet = air.wet_bulb_temperature()
    water_in = WaterFlow(temp=Q_(35, u.degC))
    solver_easy = MerkelSolver(1.0, 0.6, air, water_in, Q_(25, u.degC))
    me_easy = solver_easy.solve_me(lg_ratio=1.0)
    delta = Q_(0.5, u.delta_degC)
    t_hard_out = t_wet + delta
    solver_hard = MerkelSolver(1.0, 0.6, air, water_in, t_hard_out)
    me_hard = solver_hard.solve_me(lg_ratio=1.0)
    assert me_hard > me_easy * 3 

def test_temperature_floating():
    C = 1.0
    n = 0.6
    lg_fixed = 1.0
    dt = Q_(10, u.delta_degC)
    air_hot = AirFlow(temp=Q_(40, u.degC), humidity=Q_(40, u.perc))
    solver_hot = MerkelSolver(1.0, 0.6, air_hot, WaterFlow(), Q_(30, u.degC))
    t_in_h, t_out_h = solver_hot.estimate_temperatures(lg_fixed, dt)
    air_cold = AirFlow(temp=Q_(20, u.degC), humidity=Q_(60, u.perc))
    solver_cold = MerkelSolver(1.0, 0.6, air_cold, WaterFlow(), Q_(20, u.degC))
    t_in_c, t_out_c = solver_cold.estimate_temperatures(lg_fixed, dt)
    assert t_out_c < t_out_h 