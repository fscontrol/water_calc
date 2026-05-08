from .common import *
from .units_descriptor import Unit, QuantityType, returns, cleans, cleans_simple
from .mixins import *
from scipy.optimize import brentq

class AirFlow(UnitMagicMixin, TemporarySetMixin):
    temperature = Unit(QuantityType.TEMPERATURE)
    pressure = Unit(QuantityType.PRESSURE)
    humidity = Unit(QuantityType.HUMIDITY)
    CPA = 1006.0
    CPW = 4180.0
    CPV = 1860.0
    EVW = 2501000.0
    
    @cleans_simple(temp=u.degC, rh=u.ratio, press=u.Pa)
    def __init__(self, temp=25.0, rh=0.4, press=101325.0, flow=1.0, **kwargs):
        self.temp, self.rh, self.press, self.flow = temp, rh, press, flow
        for k, v in kwargs.items():
            setattr(self, k, v)
    
    @returns(u.degC)
    @cleans(temp=u.degC, rh=u.ratio, press=u.Pa)
    def wet_bulb_temperature(self, temp=None, rh=None, press=None):
        return psychrolib.GetTWetBulbFromRelHum(temp, rh, press)
    
    @cleans(temp=u.degC, rh=u.ratio, press=u.Pa)
    def omega(self, temp=None, rh=None, press=None):
        return  psychrolib.GetHumRatioFromRelHum(temp, rh, press)
    
    @returns(u.joule/u.kg)
    @cleans(temp=u.degC)
    def dry_air_enthalpy(self, temp=None):
        return psychrolib.GetDryAirEnthalpy(temp)

    @returns(u.joule/u.kg)
    @cleans(temp=u.degC, rh=u.ratio, press=u.Pa)
    def wet_air_enthalpy(self, temp=None, rh=None, press=None):
        omega = self.omega(temp, rh, press)
        return psychrolib.GetMoistAirEnthalpy(temp, omega)

    @returns(u.joule/u.kg)
    @cleans(temp=u.degC)
    def saturated_air_enthalpy(self, temp=None):
        omega = self.saturated_omega(temp)
        return psychrolib.GetMoistAirEnthalpy(temp, omega)

    @returns(u.g/u.kg)
    @cleans(temp=u.degC, press=u.Pa)
    def saturated_omega(self, temp=None, press=None):
        return psychrolib.GetHumRatioFromRelHum(temp, 1.0, press)
    
    @cleans(temp=u.degC, press=u.Pa)
    def lewis_factor(self, temp=None, press=None):
        w_s = self.saturated_omega(temp, press)
        return 0.865 ** (2/3) * ((1 + w_s) / (1 + 1.608 * w_s))
    
    @returns(u.joule/u.kg)
    @cleans(temp=u.degC, press=u.Pa)
    def vapor_enthalpy(self, temp=None):
        hv = self.EVW + self.CPV * temp
        return hv 
    
    @returns(u.degC)
    @cleans_simple(ha=u.joule/u.kg, tw=u.degC)
    def t_air_from_enthalpy_saturated(self, ha: float, tw: float) -> float:
        def residual(t_air):
            w_sat = psychrolib.GetHumRatioFromRelHum(t_air, 1.0, self.press)
            h_sat = psychrolib.GetMoistAirEnthalpy(t_air, w_sat)
            return h_sat - ha
        h_lo = residual(0.1)
        h_hi = residual(tw - 0.01)
        if h_lo * h_hi > 0:
            return self.temperature_from_h_w(ha, 
                psychrolib.GetHumRatioFromRelHum(
                    self.temperature_from_h_w(ha, 
                        psychrolib.GetHumRatioFromRelHum(20.0, 1.0, self.press)
                    ), 1.0, self.press
                )
            )
        return brentq(residual, 0.1, tw - 0.01, xtol=1e-4)

    @returns(u.joule/u.kg/u.degC)
    @cleans(temp=u.degC, press=u.Pa)
    def cp_saturated_air(self, temp=None, press=None):
        dt = 0.05
        w1 = psychrolib.GetHumRatioFromRelHum(temp + dt, 1.0, press)
        w2 = psychrolib.GetHumRatioFromRelHum(temp - dt, 1.0, press)
        h1 = psychrolib.GetMoistAirEnthalpy(temp + dt, w1)
        h2 = psychrolib.GetMoistAirEnthalpy(temp - dt, w2)
        return (h1 - h2) / (2 * dt)

    
    @returns(u.g/u.kg/u.degC)
    @cleans(temp=u.degC, press=u.Pa)
    def dw_sat_dT(self, temp=None, press=None):
        dt = 0.05
        w1 = psychrolib.GetHumRatioFromRelHum(temp + dt, 1.0, press)
        w2 = psychrolib.GetHumRatioFromRelHum(temp - dt, 1.0, press)
        return (w1 - w2) / (2 * dt)

    @returns(u.degC)
    @cleans(h_j_kg=u.joule/u.kg, omega_kg_kg=u.kg/u.kg)
    def temperature_from_h_w(self, h_j_kg:float, omega_kg_kg:float) -> float:
        return (h_j_kg - self.EVW * omega_kg_kg) / (self.CPA + self.CPV * omega_kg_kg)


    @returns(u.ratio)
    @cleans(temp=u.degC, omega=u.kg/u.kg, press=u.Pa)
    def rh_from_temperature_omega(self, temp: float, omega: float, press=None) -> float:
        return psychrolib.GetRelHumFromHumRatio(temp, omega, self.press)
    
    def __str__(self):
        return f"AirFlow(temp={self.temp}, rh={self.rh}, press={self.press})"
    
    @returns(u.kg/u.m**3)
    @cleans(temp=u.degC, rh=u.ratio, press=u.Pa)
    def density(self, temp=None, rh=None, press=None):
        w = self.omega(temp, rh, press)
        return psychrolib.GetMoistAirDensity(temp, w, press)
    
    @returns(u.kg)
    @cleans(temp=u.degC, rh=u.ratio, press=u.Pa)
    @cleans_simple(volume=u.m**3)
    def volume_to_dry_mass(self, volume, temp=None, rh=None, press=None):
        w = self.omega(temp, rh, press)
        rho = psychrolib.GetMoistAirDensity(temp, w, press)
        return volume * rho / (1.0 + w)
