from contextlib import contextmanager
from .common import *

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