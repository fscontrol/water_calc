from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from .common import *
from .air import AirFlow
from .water import WaterFlow
import pandas as pd
from scipy.integrate import trapezoid
from .mixins import SolverMixin, TemporarySetMixin
from .units_descriptor import cleans_simple


class PoppeSolver(SolverMixin, TemporarySetMixin):
    def poppe_system(self, tw, y):
        ha, omega = y
        h_sw = self.air_in.saturated_air_enthalpy(temp=tw)
        w_sw = self.air_in.saturated_omega(temp=tw)
        le = self.air_in.lewis_factor(temp=tw)
        hv = self.air_in.vapor_enthalpy(temp=tw)
        cp_w = self.water_in.specific_heat(temp=tw)
        d_ha_dtw = self.lg_ratio * cp_w
        t_air_dry = self.air_in.temperature_from_h_w(ha, omega)
        w_sw_air = self.air_in.omega(temp=t_air_dry, rh=1.0)
        if omega >= w_sw_air:
            t_air_sat = self.air_in.t_air_from_enthalpy_saturated(ha, tw)
            cp_as = self.air_in.cp_saturated_air(temp=t_air_sat)
            dT_air_dT_w = d_ha_dtw / cp_as
            dw_dT_air = self.air_in.dw_sat_dT(temp=t_air_sat)
            d_omega_dtw = dw_dT_air * dT_air_dT_w
        else:
            denom = (h_sw - ha) + (le - 1.0) * ((h_sw - ha) - (w_sw - omega) * hv)
            d_omega_dtw = (cp_w * self.lg_ratio * (w_sw - omega)) / max(denom, 1e-7)
        return [d_ha_dtw, d_omega_dtw]

    def solve(self):
        ha_in = self.air_in.wet_air_enthalpy()
        omega_in = self.air_in.omega()
        y0 = [ha_in, omega_in]

        try:
            sol = solve_ivp(
                self.poppe_system,
                (self.water_out.temp, self.water_in.temp),
                y0,
                t_eval=np.linspace(self.water_out.temp, self.water_in.temp, 50),
                method='Radau',
                rtol=1e-4,
                atol=1e-7,
                max_step=0.5
            )
        except ValueError as e:
            raise ValueError(f"Poppe solver failed: {e}") from None

        if sol.status != 0:
            raise ValueError(
                f"Poppe solver did not converge: {sol.message}\n"
                f"Reached tw={sol.t[-1]:.2f}°C of {self.water_in.temp:.2f}°C"
            )
        air_temperatures, air_humidities = self.decode_results(h_array=sol.y[0], omega_array=sol.y[1], press=self.air_in.press)
        zones = []
        fog_water = []
        for i, (h, w, t_a) in enumerate(zip(sol.y[0], sol.y[1], air_temperatures)):
            w_sw_air = self.air_in.omega(temp=t_a, rh=1.0)
            if w > w_sw_air:
                zones.append("fog")
                fog_water.append((w - w_sw_air))  
            else:
                zones.append("unsaturated")
                fog_water.append(0.0)

        self.profiles = pd.DataFrame({
            "water_temp_c":      sol.t,
            "air_temp_c":        air_temperatures,
            "air_rh_perc":       air_humidities*100.0,
            "air_omega_kg_kg":   sol.y[1],
            "air_enthalpy_j_kg": sol.y[0],
            "fog_kg_kg":          fog_water,   # туман в г воды на кг сухого воздуха
            "zone":              zones
        })
        omega_out = sol.y[1][-1]
        omega_in  = sol.y[1][0]

        t_air_out = air_temperatures[-1]
        w_sw_air_out = self.air_in.omega(temp=t_air_out, rh=1.0)

        if omega_out > w_sw_air_out:            # Часть omega — это туман (капли), часть — пар
            omega_vapor_out = w_sw_air_out          # пар
            omega_fog_out   = omega_out - w_sw_air_out  # туман (капли)
        else:
            omega_vapor_out = omega_out
            omega_fog_out   = 0.0

        self.evaporation  = (omega_vapor_out - omega_in) / self.lg_ratio  
        self.fog_carryover = omega_fog_out / self.lg_ratio                 
        self.total_loss   = self.evaporation + self.fog_carryover  
        self.fog_force = (self.profiles['zone'] == "fog").sum() / self.profiles['zone'].size*100.0            

        h_sw_arr = [self.air_in.saturated_air_enthalpy(temp=tw) 
                    for tw in self.profiles['water_temp_c']]
        cp_w_arr = [self.water_in.specific_heat(temp=tw) 
                    for tw in self.profiles['water_temp_c']]
        integrand = [cp / max(h_sw - ha, 1.0) 
                    for cp, h_sw, ha in zip(cp_w_arr, h_sw_arr, self.profiles['air_enthalpy_j_kg'])]
        me = trapezoid(integrand, self.profiles['water_temp_c'])       
        return float(me)
    
    @cleans_simple(press=u.Pa)
    def decode_results(self, h_array: np.array, omega_array: np.array, press: float) -> tuple[np.array, np.array]:
        t_air_list = []
        rh_list = []
        for h, w in zip(h_array, omega_array):
            t_a = self.air_in.temperature_from_h_w(h_j_kg=h, omega_kg_kg=w)
            rh = self.air_in.rh_from_temperature_omega(temp=t_a, omega=w, press=press)
            t_air_list.append(t_a)
            rh_list.append(rh)
        return np.array(t_air_list), np.array(rh_list) 