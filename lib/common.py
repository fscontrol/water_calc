import pint
import psychrolib

u = pint.UnitRegistry()
u.define('ratio = [humidity] = fraction')
u.define('perc = 0.01 * ratio = rh')
psychrolib.SetUnitSystem(psychrolib.SI)
Q_ = u.Quantity