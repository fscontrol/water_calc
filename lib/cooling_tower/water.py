from .common import *

class WaterFlow:
    def __init__(self, temp=Q_(35, u.degC), salinity=Q_(1000, u.mg/u.kg)):
        self.temp = temp
        self.salinity = salinity 

    def density(self):
        t_c = self.temp.to(u.degC).magnitude
        s = self.salinity.to(u.g/u.kg).magnitude
        rho_w = (999.842594 + 6.793952e-2 * t_c - 9.095290e-3 * t_c**2 + 
                 1.001685e-4 * t_c**3 - 1.120083e-6 * t_c**4 + 6.536332e-9 * t_c**5)
        a = 8.24493e-1 - 4.0899e-3 * t_c + 7.6438e-5 * t_c**2 - 8.2467e-7 * t_c**3 + 5.3875e-9 * t_c**4
        b = -5.72466e-3 + 1.0227e-4 * t_c - 1.6546e-6 * t_c**2
        c = 4.8314e-4
        rho_brine = rho_w + a*s + b*(s**1.5) + c*(s**2)
        return Q_(rho_brine, u.kg / u.m**3)

    def viscosity(self):
        t_c = self.temp.to(u.degC).magnitude
        s = self.salinity.to(u.g/u.kg).magnitude
        mu_w = 0.001 * (1 + 0.03368 * t_c + 0.000221 * t_c**2)**-1
        a_coeff = 0.005 
        b_coeff = 0.145
        s_norm = s / 58.44    
        mu_brine = mu_w * (1 + a_coeff * (s_norm**0.5) + b_coeff * s_norm)
        return Q_(mu_brine, u.Pa * u.s)

    def latent_heat(self):
        t_c = self.temp.to(u.degC).magnitude
        s = self.salinity.to(u.g/u.kg).magnitude
        lv_w = 2501000 - 2360 * t_c 
        s_fraction = s / 1000.0
        lv_brine = lv_w * (1 - s_fraction)
        return Q_(lv_brine, u.J / u.kg)
    def specific_heat(self, temperature = None, salinity=None):
        if temperature is None:
            t_c = self.temp.to(u.degC).magnitude
        else:
            t_c = temperature.to(u.degC).magnitude
        if salinity is None:
            s = self.salinity.to(u.g/u.kg).magnitude
        else:
            s = salinity.to(u.g/u.kg).magnitude
        cp_fresh = 4185.5 + 0.1 * t_c
        cp_brine = cp_fresh - s * (5.2 - 0.02 * t_c)
        return Q_(cp_brine, u.J / (u.kg * u.delta_degC))