from enum import Enum, auto
import pint 
from functools import wraps
import inspect
import pint
u = pint.UnitRegistry()
Q_ = u.Quantity

u = pint.UnitRegistry()
u.define('ratio = [humidity] = fraction')
u.define('perc = 0.01 * ratio = rh')
Q_ = u.Quantity

class QuantityType(Enum):
    """Типы физических величин"""
    TEMPERATURE = auto()
    PRESSURE = auto()
    HUMIDITY = auto()
    ENTHALPY = auto()
    SPECIFIC_HEAT = auto()
    SPECIFIC_ENTHALPY = auto()
    FLOW = auto()
    TDS = auto()

def quantity_value(a, b):
    return dict(base_unit=a, field=b)

class Unit:
    CONFIG = {
        QuantityType.TEMPERATURE: quantity_value(u.degC, "temp"),
        QuantityType.PRESSURE: quantity_value(u.Pa, "press"),
        QuantityType.HUMIDITY: quantity_value(u.ratio,"rh"),
        QuantityType.ENTHALPY: quantity_value(u.joule, "H"),
        QuantityType.SPECIFIC_HEAT: quantity_value(u.joule/u.kg/u.degC, "cp"),
        QuantityType.SPECIFIC_ENTHALPY: quantity_value(u.joule/u.kg, "h"),
        QuantityType.TDS: quantity_value(u.g/u.kg, "salinity"),
    }
    def __init__(self, quantity_type: QuantityType):
        self.quantity_type = quantity_type
        self.config = self.CONFIG[quantity_type]
    
    def __get__(self, instance, owner):
        if instance is None:
            return self       
        value = getattr(instance, self.config["field"])
        return Q_(value, self.config['base_unit'])
    
    def __set__(self, instance, value):
        if hasattr(value, 'magnitude'):
            val = value.to(self.config['base_unit']).magnitude
        else:
            val = float(value)
        setattr(instance, self.config["field"], val)

def returns(unit):
    def decorator(func):
        func._return_unit = unit
        return func
    return decorator

def cleans(**units_map):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            self = args[0]
            sig = inspect.signature(func)
            bound = sig.bind(*args, **kwargs)
            bound.apply_defaults()
            for arg_name, arg_value in bound.arguments.items():
                if arg_name == 'self':
                    continue
                if arg_name in units_map:
                    bound.arguments[arg_name] = clean_param(self, arg_name, arg_value, units_map[arg_name])
            return func(*bound.args, **bound.kwargs)
        return wrapper
    return decorator

def clean_param(self, param_name, param_value, unit):
    if param_value is None:
        if hasattr(self, param_name):
            return getattr(self, param_name)
        else:
            raise Exception("no param in self, no param given")
    return clean_param_simple(self, param_name, param_value, unit)

def cleans_simple(**units_map):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            self = args[0]
            sig = inspect.signature(func)
            bound = sig.bind(*args, **kwargs)
            bound.apply_defaults()
            for arg_name, arg_value in bound.arguments.items():
                if arg_name == 'self':
                    continue
                if arg_name in units_map:
                    bound.arguments[arg_name] = clean_param_simple(self, arg_name, arg_value, units_map[arg_name])
            return func(*bound.args, **bound.kwargs)
        return wrapper
    return decorator

def clean_param_simple(self, param_name, param_value, unit):
    try:
        return param_value.to(unit).magnitude
    except:
        return float(param_value)