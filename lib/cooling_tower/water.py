from .common import *
from .units_descriptor import QuantityType, Unit, cleans, returns, cleans_simple
from contextlib import contextmanager
from .mixins import UnitMagicMixin, TemporarySetMixin

class WaterFlow(UnitMagicMixin, TemporarySetMixin):
    temperature = Unit(QuantityType.TEMPERATURE)
    tds = Unit(QuantityType.TDS)
    
    @cleans_simple(temp=u.degC, salinity=u.g/u.kg)
    def __init__(self, temp=35.0, salinity=1.0, **kwargs):
        self.temp, self.salinity = temp, salinity 
        for k, v in kwargs.items():
            setattr(self, k, v)
    
    @returns(u.kg/u.m**3)
    @cleans(temp=u.degC, salinity=u.g/u.kg)
    def density(self, temp=None, salinity=None):
        rho_w = (999.842594 + 6.793952e-2 * temp - 9.095290e-3 * temp**2 + 
                 1.001685e-4 * temp**3 - 1.120083e-6 * temp**4 + 6.536332e-9 * temp**5)
        a = 8.24493e-1 - 4.0899e-3 * temp + 7.6438e-5 * temp**2 - 8.2467e-7 * temp**3 + 5.3875e-9 * temp**4
        b = -5.72466e-3 + 1.0227e-4 * temp - 1.6546e-6 * temp**2
        c = 4.8314e-4
        rho_brine = rho_w + a*salinity + b*(salinity**1.5) + c*(salinity**2)
        return rho_brine

    @returns(u.Pa*u.s)
    @cleans(temp=u.degC, salinity=u.g/u.kg)
    def viscosity(self, temp=None, salinity=None):
        t_k = temp + 273.15
        mu_w = np.exp(5.495921e5 / t_k**2 - 1.66779e3 / t_k - 7.612821)
        A = 1.541 + 1.998e-2 * temp + 9.585e-5 * temp**2
        B = 7.974 - 7.561e-2 * temp + 4.724e-4 * temp**2
        C = 1.0 / (0.001 * (temp + 91.012) + 0.003380 * temp**2)
        mu_rel = 1 + A * salinity + B * salinity**2 + C * salinity**3    
        mu = mu_w * mu_rel
        return mu

    @returns(u.joule/u.kg)
    @cleans(temp=u.degC, salinity=u.g/u.kg)
    def latent_heat(self, temp=None, salinity=None):
        lv_w = 2501000 - 2360 * temp 
        s_fraction = salinity / 1000.0
        lv_brine = lv_w * (1 - s_fraction)
        return lv_brine
    
    @returns(u.joule/u.kg/u.degC)
    @cleans(temp=u.degC, salinity=u.g/u.kg)
    def specific_heat(self, temp = None, salinity=None):
        cp_fresh = 4185.5 + 0.1 * temp
        cp_brine = cp_fresh - salinity * (5.2 - 0.02 * temp)
        return cp_brine