import pytest
import numpy as np
from ..lib import *

def test_poppe_solver_execution():
    """Check if the solver runs and returns a non-empty DataFrame."""
    air_in = AirFlow(temp=Q_(25, u.degC), humidity=Q_(60, u.perc))
    water_in = WaterFlow(temp=Q_(40, u.degC))
    t_out = Q_(30, u.degC)
    lg = 1.0
    
    solver = PoppeSolver(air_in, water_in, t_out, lg)
    df = solver.solve()
    
    assert not df.empty
    assert 'air_enthalpy_j_kg' in df.columns
    assert 'water_temp_c' in df.columns

def test_poppe_physics_trends():
    """Air enthalpy should increase as it moves through the tower."""
    air_in = AirFlow(temp=Q_(20, u.degC), humidity=Q_(50, u.perc))
    water_in = WaterFlow(temp=Q_(45, u.degC))
    t_out = Q_(32, u.degC)
    
    solver = PoppeSolver(air_in, water_in, t_out, lg_ratio=1.0)
    df = solver.solve()
    
    h_start = df['air_enthalpy_j_kg'].iloc[0]
    h_end = df['air_enthalpy_j_kg'].iloc[-1]
    
    # In a cooling tower, air absorbs heat
    assert h_end > h_start

def test_poppe_mass_balance():
    """Evaporation should be positive in standard conditions."""
    air_in = AirFlow(temp=Q_(30, u.degC), humidity=Q_(30, u.perc))
    water_in = WaterFlow(temp=Q_(45, u.degC))
    solver = PoppeSolver(air_in, water_in, Q_(35, u.degC), lg_ratio=1.1)
    df = solver.solve()
    
    # Check attributes stored in metadata
    assert df.attrs['evaporation_kg_kg'] > 0
    assert df.attrs['total_water_loss_kg_kg'] >= df.attrs['evaporation_kg_kg']
