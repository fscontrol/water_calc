import pint
import psychrolib
import numpy as np

u = pint.UnitRegistry()
u.define('ratio = [humidity] = fraction')
u.define('perc = 0.01 * ratio = rh')
Q_ = u.Quantity

class Temperature:
    def __get__(self, instance, owner):
        if instance is None:
            return self
        return Q_(instance.temp, u.degC)
    def __set__(self, instance, value):
        if hasattr(value, 'magnitude'):
            instance.temp = value.to(u.degC).magnitude
        else:
            instance.temp = float(value)

class Air:
    temperature = Temperature() 
    def __init__(self, t = 25.0):
        self.temp = t

air = Air()
air.temperature = Q_(44, u.degC)
print(air.temp)
