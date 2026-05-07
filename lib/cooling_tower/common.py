import pint
import psychrolib
import numpy as np

u = pint.UnitRegistry()
u.define('ratio = [humidity] = fraction')
u.define('perc = 0.01 * ratio = rh')
psychrolib.SetUnitSystem(psychrolib.SI)
Q_ = u.Quantity

class TemperatureUnit:
    units = dict(
        celsius = [1, 0],
        farenheit = [1.8, 32],
        kelvin = [1, 273.15]
    )
    def __set_name__(self, name):
        self.name = name
    def __set__(self, instance, value):
        instance.temp = (value - self.units[self.name][1])/self.units[self.name][0]
    
    def __get__(self, instance, owner):
        return instance.temp*self.units[self.name][0] + self.units[self.name][1]
