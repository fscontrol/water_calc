import numpy as np
from scipy.integrate import quad
from .air import AirFlow
from .water import WaterFlow
from .common import u, Q_
from scipy.optimize import fsolve

class MerkelSolver:
    def __init__(self, C, n, air_in: AirFlow, water_in: WaterFlow, water_temp_out: Q_):
        self.air_in = air_in
        self.water_in = water_in
        self.t_out = water_temp_out
        self.C = C
        self.n = n

    def get_air_enthalpy_at_tw(self, t_w_current, lg_ratio):
        h_a_in = self.air_in.wet_air_enthalpy()
        cp_w = self.water_in.specific_heat()
        delta_t = t_w_current - self.t_out.to(u.degC).magnitude
        h_a = h_a_in.magnitude + lg_ratio * cp_w.magnitude * delta_t
        return h_a

    def solve_me(self, lg_ratio):
        cp_w = self.water_in.specific_heat().magnitude
        t_in_c = self.water_in.temp.to(u.degC).magnitude
        t_out_c = self.t_out.to(u.degC).magnitude
        def integrand(tw):
            cp_w = self.water_in.specific_heat(temperature=Q_(tw, u.degC)).magnitude
            h_sw = self.air_in.saturated_air_enthalpy(Q_(tw, u.degC)).magnitude
            h_a = self.get_air_enthalpy_at_tw(tw, lg_ratio)
            diff = h_sw - h_a
            return cp_w / max(diff, 1e-4)
        me_number, error = quad(integrand, t_out_c, t_in_c)
        return me_number
    
    def target_me(self, lg_ratio):
        return self.C * (lg_ratio**(-self.n))

    def find_operating_lg(self):
        def objective(lg_ratio):
            me_req = self.solve_me(lg_ratio)
            me_avail = self.target_me(lg_ratio)
            return me_req - me_avail
        lg_final = fsolve(objective, 1.0)
        return lg_final[0]
    
    def estimate_temperatures(self,lg_ratio, delta_t=Q_(10, u.delta_degC)):
        dt_val = delta_t.to(u.delta_degC).magnitude
        t_wet_q = self.air_in.wet_bulb_temperature()
        t_wet = t_wet_q.magnitude
        def objective(t_out_guess):
            current_t_out = float(t_out_guess[0])
            if current_t_out <= t_wet:
                return 1e6 * (t_wet - current_t_out + 1)
            current_t_in = current_t_out + dt_val          
            temp_solver = MerkelSolver(self.C, self.n,
                self.air_in, 
                WaterFlow(temp=Q_(current_t_in, u.degC)), 
                Q_(current_t_out, u.degC)
            )       
            try:
                val = temp_solver.solve_me(lg_ratio) - self.target_me(lg_ratio)
                return val
            except Exception:
                return 1e6

        t_out_final = fsolve(objective, t_wet + 3.0, xtol=1e-4)
        t_out_res = Q_(float(t_out_final[0]), u.degC)
        t_in_res = t_out_res + delta_t
        return t_in_res, t_out_res