import pint
import psychrolib
import numpy as np
#from .units_register import u, Q_

psychrolib.SetUnitSystem(psychrolib.SI)
import pint

u = pint.UnitRegistry()
Q_ = u.Quantity

u = pint.UnitRegistry()
u.define('ratio = [humidity] = fraction')
u.define('perc = 0.01 * ratio = rh')
Q_ = u.Quantity