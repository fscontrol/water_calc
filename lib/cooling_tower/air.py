from .common import *
class AirFlow:
    def __init__(self, temp=Q_(25, u.degC), humidity=Q_(40.0, u.perc), pressure=Q_(1, u.atm)):
        self.temp = temp.to(u.degC).magnitude
        self.humidity = humidity.to(u.ratio).magnitude
        self.pressure = pressure.to(u.Pa).magnitude
    def wet_bulb_temperature(self):
        t_wet = psychrolib.GetTWetBulbFromRelHum(self.temp, self.humidity, self.pressure)
        return Q_(t_wet, u.degC)
    def omega(self):
        w = psychrolib.GetHumRatioFromRelHum(self.temp, self.humidity, self.pressure)
        return Q_(w, u.ratio)
    
    def dry_air_enthalpy(self):
        h_dry = psychrolib.GetDryAirEnthalpy(self.temp)
        return Q_(h_dry, u.J / u.kg)
    def wet_air_enthalpy(self):
        w = self.omega().magnitude
        h_wet = psychrolib.GetMoistAirEnthalpy(self.temp, w)
        return Q_(h_wet, u.J / u.kg)

    def saturated_air_enthalpy(self, t_water):
            tw_c = t_water.to(u.degC).magnitude
            w_sat = psychrolib.GetHumRatioFromRelHum(tw_c, 1.0, self.pressure)
            h_sat = psychrolib.GetMoistAirEnthalpy(tw_c, w_sat)
            return Q_(h_sat, u.J / u.kg)

    def saturated_omega(self, t_water):
        tw_c = t_water.to(u.degC).magnitude
        w_sat = psychrolib.GetHumRatioFromRelHum(tw_c, 1.0, self.pressure)
        return Q_(w_sat, u.ratio)
    
    def lewis_factor(self):
        w_s = self.saturated_omega(Q_(self.temp, u.degC)).magnitude
        le = 0.865 ** (2/3) * ((1 + w_s) / (1 + 1.608 * w_s))
        return le
    
    