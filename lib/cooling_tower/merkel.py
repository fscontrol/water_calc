from scipy.integrate import quad
from scipy.optimize import fsolve, brentq # Добавили brentq
from .air import AirFlow
from .water import WaterFlow
from .common import *
class MerkelSolver:
    def __init__(self, C, n, air_in: AirFlow, water_in: WaterFlow, water_temp_out: Q_):
        self.air_in = air_in
        self.water_in = water_in
        self.t_out = water_temp_out
        self.C = C
        self.n = n

    def get_air_enthalpy_at_tw(self, t_w_current, lg_ratio):
        h_a_in = self.air_in.wet_air_enthalpy()
        cp_w = self.water_in.specific_heat(Q_(t_w_current, u.degC)).magnitude
        delta_t = t_w_current - self.t_out.to(u.degC).magnitude
        h_a = h_a_in.magnitude + lg_ratio * cp_w * delta_t
        return h_a

    def solve_me(self, lg_ratio):
        twb = self.air_in.wet_bulb_temperature()
        if twb > self.t_out:
            raise Exception(f"t_out has to be upper then wet bulb temperature, Twb = {twb}")
        t_in_c  = self.water_in.temp.to(u.degC).magnitude
        t_out_c = self.t_out.to(u.degC).magnitude
        
        def integrand(tw):
            cp_w = self.water_in.specific_heat(
                temperature=Q_(tw, u.degC)).magnitude
            h_sw = self.air_in.saturated_air_enthalpy(Q_(tw, u.degC)).magnitude
            h_a  = self.get_air_enthalpy_at_tw(tw, lg_ratio)
            diff = h_sw - h_a
            return cp_w / max(diff, 0.1)
        
        me_number, _ = quad(integrand, t_out_c, t_in_c, limit=100)
        return me_number
    
    def target_me(self, lg_ratio):
        return self.C * (lg_ratio**(-self.n))

    def find_operating_lg(self):
        def objective(lg_ratio_arr):
            lg = lg_ratio_arr[0]
            me_req = self.solve_me(lg)
            me_avail = self.target_me(lg)
            return me_req - me_avail
            
        lg_final = fsolve(objective, 1.0, xtol=1e-4)
        return float(lg_final[0])
    
    def estimate_temperatures(self, lg_ratio, delta_t=Q_(10, u.delta_degC)):
        dt_val = delta_t.to(u.delta_degC).magnitude
        twb = self.air_in.wet_bulb_temperature().magnitude
        me_avail = self.target_me(lg_ratio)
        error_message = None

        def objective(t_out_guess):
            current_t_out = float(t_out_guess[0] if hasattr(t_out_guess, "__len__") else t_out_guess)
            current_t_in = current_t_out + dt_val          
            temp_solver = MerkelSolver(
                self.C, self.n, self.air_in, 
                WaterFlow(temp=Q_(current_t_in, u.degC)), 
                Q_(current_t_out, u.degC)
            )       
            return temp_solver.solve_me(lg_ratio) - me_avail
        a = twb + 2.0
        b = 90.0 - dt_val  
        fa = objective(a)
        fb = objective(b)           
        if fa * fb < 0:
            t_out_final = brentq(objective, a, b, xtol=1e-4)
        else:
            if fa < 0:
                error_message = f"Coolign tower Merkel is higher than required Tlow={a}, LG = {lg_ratio} ME = {me_avail}"
                t_out = Q_(a,u.degC)
                return t_out + delta_t, t_out, error_message
            else:
                error_message = f"Cooling tower Merkel is lower than required Tlow={a}, LG = {lg_ratio} ME = {me_avail}"
                raise Exception(error_message)

        t_out_res = Q_(float(t_out_final), u.degC)
        t_in_res = t_out_res + delta_t
        return t_in_res, t_out_res, error_message