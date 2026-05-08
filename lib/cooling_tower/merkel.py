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


    