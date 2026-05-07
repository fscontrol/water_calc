from scipy.integrate import quad
from scipy.optimize import fsolve, brentq # Добавили brentq
from .air import AirFlow
from .water import WaterFlow
from .common import *
from contextlib import contextmanager
from .mixins import TemporarySetMixin, UnitMagicMixin, SolverMixin
from .units_descriptor import Unit, cleans_simple, returns

class MerkelSolver(SolverMixin, TemporarySetMixin):
    @returns(u.joule/u.kg)
    @cleans_simple(t_w_current=u.degC)
    def air_enthalpy_at_tw(self, t_w_current, lg_ratio=None):
        lg = lg_ratio if lg_ratio is not None else self.lg_ratio
        h_a_in = self.air_in.wet_air_enthalpy()
        cp_w = self.water_in.specific_heat(t_w_current)
        delta_t = t_w_current - self.water_out.temp
        h_a = h_a_in + lg * cp_w * delta_t
        return h_a

    def solve(self, lg_ratio=None):
        lg = lg_ratio if lg_ratio is not None else self.lg_ratio
        twb = self.air_in.wet_bulb_temperature()
        if twb > self.water_out.temp:
            raise Exception(f"t_out has to be upper then wet bulb temperature, Twb = {twb}")
        t_in_c  = self.water_in.temp
        t_out_c = self.water_out.temp
        
        def integrand(tw):
            cp_w = self.water_in.specific_heat(tw)
            h_sw = self.air_in.saturated_air_enthalpy(tw)
            h_a  = self.air_enthalpy_at_tw(tw, lg)
            diff = h_sw - h_a
            return cp_w / max(diff, 0.1)
        
        me_number, _ = quad(integrand, t_out_c, t_in_c, limit=100)
        return me_number

    def find_operating_lg(self):
        def objective(lg_ratio_arr):
            lg = lg_ratio_arr[0]
            me_req = self.solve(lg)
            me_avail = self.target_me(lg)
            return me_req - me_avail
            
        lg_final = fsolve(objective, 1.0, xtol=1e-4)
        return float(lg_final[0])
    
    @cleans_simple(delta_t=u.delta_degC)
    def estimate_temperatures(self, lg_ratio=None, delta_t=10.0, return_units=False):
        lg = lg_ratio if lg_ratio is not None else self.lg_ratio
        dt_val = delta_t
        twb = self.air_in.wet_bulb_temperature()
        me_avail = self.target_me(lg)
        error_message = None

        def objective(t_out_guess):
            current_t_out = float(t_out_guess[0] if hasattr(t_out_guess, "__len__") else t_out_guess)
            current_t_in = current_t_out + dt_val          
            temp_solver = MerkelSolver(
                C=self.C, n=self.n, air_in=self.air_in, 
                water_in=WaterFlow(temp=current_t_in), 
                water_out=WaterFlow(temp=current_t_out),
                lg_ratio=lg
            )       
            return temp_solver.solve() - me_avail
        a = twb + 2.0
        b = 90.0 - dt_val  
        fa = objective(a)
        fb = objective(b)           
        if fa * fb < 0:
            t_out_final = brentq(objective, a, b, xtol=1e-4)
        else:
            if fa < 0:
                error_message = f"Coolign tower Merkel is higher than required Tlow={a}, LG = {lg_ratio} ME = {me_avail}"
                t_out = a
                return t_out + delta_t, t_out, error_message
            else:
                error_message = f"Cooling tower Merkel is lower than required Tlow={a}, LG = {lg_ratio} ME = {me_avail}"
                raise Exception(error_message)

        t_out_res = t_out_final
        t_in_res = t_out_res + delta_t
        if return_units:
            return Q_(t_in_res, u.degC), Q_(t_out_res, u.degC), error_message
        else:
            return t_in_res, t_out_res, error_message