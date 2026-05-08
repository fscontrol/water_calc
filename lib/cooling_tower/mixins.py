from __future__ import annotations
from contextlib import contextmanager
from .common import *
from .units_descriptor import cleans_simple
from scipy.optimize import brentq


class UnitMagicMixin:
    def __getattr__(self, name):
        if name.startswith("calc_"):
            cache_attr = f"_{name}_cached"
            if hasattr(self, cache_attr):
                return getattr(self, cache_attr)
            method_name = name[5:]
            method = getattr(self, method_name, None)
            if method is None or not hasattr(method, '_return_unit'):
                raise AttributeError(f"'{method_name}' not found or has no unit info")
            def wrapper(unit=None, *args, **kwargs):
                result = method(*args, **kwargs)
                q = Q_(result, method._return_unit)
                return q.to(unit) if unit is not None else q
            setattr(self, cache_attr, wrapper)
            return wrapper
        raise AttributeError(name)

class TemporarySetMixin:
    @contextmanager
    def temporary_set(self, **kwargs):
        old_values = {}
        for key, value in kwargs.items():
            if hasattr(self, key):
                old_values[key] = getattr(self, key)
                setattr(self, key, value)
            else:
                raise AttributeError(f"{self.__class__.__name__} has no attribute '{key}'")
        try:
            yield self
        finally:
            for key, value in old_values.items():
                setattr(self, key, value)

class SolverMixin:
    def __init__(self, air_in: AirFlow, water_in: WaterFlow, water_out: WaterFlow, lg_ratio: float=1.0, C=1.0, n=0.6):
        self.air_in = air_in
        self.water_in = water_in
        self.water_out = water_out
        self.lg_ratio = lg_ratio
        self.C = C
        self.n = n

    def _validate_temperatures(self):
        twb = self.air_in.wet_bulb_temperature()
        if twb > self.water_out.temp:
            raise ValueError(f"T_out must be above wet bulb: Twb={twb:.1f}°C")
    
    def _water_at(self, temp_c):
        """Создать WaterFlow с текущей солёностью при заданной температуре"""
        return WaterFlow(temp=temp_c, salinity=self.water_in.salinity)
    
    def target_me(self, lg_ratio=None):
        lg = lg_ratio if lg_ratio is not None else self.lg_ratio
        return self.C * (lg**(-self.n))

    @cleans_simple(delta_t=u.delta_degC)
    def estimate_temperatures(self, lg_ratio=None, delta_t=10.0, target_me=None, return_units=False):
        lg = lg_ratio if lg_ratio is not None else self.lg_ratio
        dt_val = delta_t
        twb = self.air_in.wet_bulb_temperature()
        me_avail = target_me if target_me is not None else self.target_me(lg)
        
        def objective(t_out_guess):
            current_t_out = float(t_out_guess[0] if hasattr(t_out_guess, "__len__") else t_out_guess)
            current_t_in = current_t_out + dt_val 
            with self.water_in.temporary_set(temp=current_t_in):
                with self.water_out.temporary_set(temp=current_t_out):   
                    temp_solver = self.__class__(
                        C=self.C, n=self.n, air_in=self.air_in, 
                        water_in=self.water_in, water_out=self.water_out, 
                        lg_ratio=lg
                    )       
                    return temp_solver.solve() - me_avail
        a = twb + 2.0
        b = 90.0 - dt_val
        flag = False
        while not flag:
            try:
                fa = objective(a)
                flag = True 
            except ValueError:
                a += 0.1
            if a > b:
                raise Exception("For this Me and LG ratio delta T is too high or T_cold is too low")
        fb = objective(b)           
        if fa * fb < 0:
            t_out_final = brentq(objective, a, b, xtol=1e-4)
        else:
            if fa < 0:
                raise Exception(f"Coolign tower Merkel is higher than required Tlow={a}, LG = {lg_ratio} ME = {me_avail}")
            else:
                raise Exception(f"Cooling tower Merkel is lower than required Tlow={a}, LG = {lg_ratio} ME = {me_avail}")

        t_out_res = t_out_final
        t_in_res = t_out_res + delta_t
        if return_units:
            return Q_(t_in_res, u.degC), Q_(t_out_res, u.degC)
        else:
            return t_in_res, t_out_res
    def find_operating_lg(self):
        def objective(lg):
            lg = float(lg)
            me_req = self.solve(lg)
            me_avail = self.target_me(lg)
            return me_req - me_avail
        a, b = 0.1, 10.0
        f_a, f_b = objective(a), objective(b)
        if f_a * f_b < 0:
            lg_final = brentq(objective, a, b, xtol=1e-4)
        else:
            raise Exception(f"No solution found for LG ratio between {a} and {b}")
        return float(lg_final)