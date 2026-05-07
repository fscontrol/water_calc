import pint
from ..lib import *
import pytest as pt

def test_target_me_power_law():
    air = AirFlow(temp=25.0, rh=0.5)
    w_in = WaterFlow(temp=40.0)
    w_out = WaterFlow(temp=30.0)
    solver = MerkelSolver(air, w_in, w_out, 1.0, 2.0, 0.6)
    assert solver.target_me(1.0) == pt.approx(2.0)
    assert solver.target_me(0.5) > solver.target_me(1.0)

def test_find_operating_lg_matches_curve():
    air = AirFlow(temp=25.0, rh=0.5)
    w_in = WaterFlow(temp=42.0)
    w_out = WaterFlow(temp=32.0)
    solver = MerkelSolver(air, w_in, w_out, 1.0, 1.2, 0.58)
    lg = solver.find_operating_lg()
    assert lg > 0
    assert solver.solve(lg) == pt.approx(solver.target_me(lg), rel=1e-3)

def test_solve_raises_when_out_below_wet_bulb():
    air = AirFlow(temp=30.0, rh=0.85)
    twb = air.wet_bulb_temperature()
    twb_c = float(getattr(twb, "magnitude", twb))
    w_in = WaterFlow(temp=45.0)
    w_out = WaterFlow(temp=twb_c - 5.0)
    solver = MerkelSolver(air, w_in, w_out, 1.0, 1.0, 0.6)
    with pt.raises(Exception, match="wet bulb"):
        solver.solve()

def test_air_enthalpy_at_outlet_equals_inlet_when_delta_t_zero():
    air = AirFlow(temp=25.0, rh=0.5)
    w_in = WaterFlow(temp=40.0)
    t_out = 30.0
    w_out = WaterFlow(temp=t_out)
    solver = MerkelSolver(air, w_in, w_out, 1.1, 1.0, 0.6)
    h_at_tout = solver.air_enthalpy_at_tw(t_out)
    h_in = air.wet_air_enthalpy()
    assert h_at_tout == pt.approx(h_in, rel=1e-5)

def test_estimate_temperatures_return_units():
    air = AirFlow(temp=25.0, rh=0.5)
    solver = MerkelSolver(air, WaterFlow(temp=35.0), WaterFlow(temp=25.0), 1.0, 1.0, 0.6)
    t_in, t_out, err = solver.estimate_temperatures(1.0, Q_(8, u.delta_degC), return_units=True)
    assert err is None
    assert hasattr(t_in, "magnitude")
    assert t_in.magnitude == pt.approx(t_out.magnitude + 8.0, abs=0.01)

def test_estimate_temperatures_delta_as_quantity_matches_scalar():
    air = AirFlow(temperature=Q_(40, u.degC), humidity=Q_(40, u.perc))
    solver = MerkelSolver(air, WaterFlow(), WaterFlow(temperature=Q_(30, u.degC)), 1.0, 1.0, 0.6)
    dt_mag = Q_(10, u.delta_degC).magnitude
    t_in1, t_out1, _ = solver.estimate_temperatures(1.0, dt_mag)
    t_in2, t_out2, _ = solver.estimate_temperatures(1.0, Q_(10, u.delta_degC))
    assert t_in1 == pt.approx(t_in2)
    assert t_out1 == pt.approx(t_out2)

def test_merkel_logic():
    air = AirFlow(temp=25.0, humidity=Q_(50, u.perc))
    water_in = WaterFlow(temperature=Q_(40, u.degC)) # Вход 40
    water_out = WaterFlow(temperature= Q_(35, u.degC))  
    solver = MerkelSolver(air, water_in, water_out, 1.0, 1.0, 0.6)
    me_val = solver.solve()
    assert me_val > 0
    water_out.temperature = Q_(30, u.degC)
    solver_colder = MerkelSolver(air, water_in, water_out, 1.0, 1.0, 0.6)
    me_colder = solver_colder.solve()
    
    assert me_colder > me_val

def test_merkel_standard_case():
    air_in = AirFlow(temperature=Q_(25, u.degC), humidity=Q_(50, u.perc))
    water_in = WaterFlow(temperature=Q_(40, u.degC))
    water_out = WaterFlow(temperature=Q_(28, u.degC))
    solver = MerkelSolver(air_in, water_in, water_out, 1.0, 1.0, 0.6)
    me_val = solver.solve()
    assert 1.0 < me_val < 1.5
    print(f"\nРассчитанное число Меркеля: {me_val:.4f}")

def test_merkel_limit_approach():
    t_dry = Q_(25, u.degC)
    rh = Q_(60, u.perc)
    air = AirFlow(temperature=t_dry, humidity=rh)
    t_wet = air.wet_bulb_temperature()
    water_in = WaterFlow(temperature=Q_(55, u.degC))
    water_out = WaterFlow(temp=45.0)
    solver_easy = MerkelSolver(air, water_in, water_out, 1.0, 1.0, 0.6)
    me_easy = solver_easy.solve()
    delta = 10.0
    t_hard_out = t_wet + delta
    water_out.temp = t_hard_out
    solver_hard = MerkelSolver(air, water_in, water_out, 1.0, 1.0, 0.6)
    me_hard = solver_hard.solve()
    assert me_hard > me_easy * 3 

def test_temperature_floating():
    lg_fixed = 1.0
    dt = Q_(10, u.delta_degC).magnitude
    air_hot = AirFlow(temperature=Q_(40, u.degC), humidity=Q_(40, u.perc))
    solver_hot = MerkelSolver(air_hot, WaterFlow(), WaterFlow(temperature=Q_(30, u.degC)), 1.0, 1.0, 0.6)
    t_in_h, t_out_h, error_message_h = solver_hot.estimate_temperatures(lg_fixed, dt)
    air_cold = AirFlow(temperature=Q_(20, u.degC), humidity=Q_(60, u.perc))
    solver_cold = MerkelSolver(air_cold, WaterFlow(), WaterFlow(temperature=Q_(20, u.degC)), 1.0, 1.0, 0.6)
    t_in_c, t_out_c, error_message_c = solver_cold.estimate_temperatures(lg_fixed, dt)
    assert t_out_c < t_out_h 
    assert error_message_h is None
    assert error_message_c is None
