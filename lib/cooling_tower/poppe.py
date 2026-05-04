from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from .common import *
from .air import AirFlow
from .water import WaterFlow
import pandas as pd
from scipy.integrate import trapezoid

class PoppeSolver:
    def __init__(self, air_in: AirFlow, water_in: WaterFlow, water_temp_out: Q_, lg_ratio: float):
        self.air_in = air_in
        self.water_in = water_in
        self.t_out = water_temp_out.to(u.degC).magnitude
        self.t_in = water_in.temp.to(u.degC).magnitude
        self.lg_ratio = lg_ratio

    def _t_air_from_ha_saturated(self, ha: float, tw: float) -> float:
        def residual(t_air):
            w_sat = psychrolib.GetHumRatioFromRelHum(t_air, 1.0, self.air_in.pressure)
            h_sat = psychrolib.GetMoistAirEnthalpy(t_air, w_sat)
            return h_sat - ha
        h_lo = residual(0.1)
        h_hi = residual(tw - 0.01)
        if h_lo * h_hi > 0:
            return AirFlow.calc_temperature_from_h_w(ha, 
                psychrolib.GetHumRatioFromRelHum(
                    AirFlow.calc_temperature_from_h_w(ha, 
                        psychrolib.GetHumRatioFromRelHum(20.0, 1.0, self.air_in.pressure)
                    ), 1.0, self.air_in.pressure
                )
            )
        return brentq(residual, 0.1, tw - 0.01, xtol=1e-4)

    def _cp_saturated_air(self, t_air: float) -> float:
        dt = 0.05
        w1 = psychrolib.GetHumRatioFromRelHum(t_air + dt, 1.0, self.air_in.pressure)
        w2 = psychrolib.GetHumRatioFromRelHum(t_air - dt, 1.0, self.air_in.pressure)
        h1 = psychrolib.GetMoistAirEnthalpy(t_air + dt, w1)
        h2 = psychrolib.GetMoistAirEnthalpy(t_air - dt, w2)
        return (h1 - h2) / (2 * dt)

    def _dw_sat_dT(self, t_air: float) -> float:
        dt = 0.05
        w1 = psychrolib.GetHumRatioFromRelHum(t_air + dt, 1.0, self.air_in.pressure)
        w2 = psychrolib.GetHumRatioFromRelHum(t_air - dt, 1.0, self.air_in.pressure)
        return (w1 - w2) / (2 * dt)

    def poppe_system(self, tw, y):
        ha, omega = y
        h_sw = self.air_in.saturated_air_enthalpy(Q_(tw, u.degC)).magnitude
        w_sw = self.air_in.saturated_omega(Q_(tw, u.degC)).magnitude
        le = self.air_in.lewis_factor(Q_(tw, u.degC))
        hv = 2501000 + 1860 * tw
        cp_w = self.water_in.specific_heat(temperature=Q_(tw, u.degC)).magnitude
        d_ha_dtw = self.lg_ratio * cp_w
        t_air_dry = AirFlow.calc_temperature_from_h_w(ha, omega)
        w_sw_air = psychrolib.GetHumRatioFromRelHum(t_air_dry, 1.0, self.air_in.pressure)
        if omega >= w_sw_air:
            t_air_sat = self._t_air_from_ha_saturated(ha, tw)
            cp_as = self._cp_saturated_air(t_air_sat)
            dT_air_dT_w = d_ha_dtw / cp_as
            dw_dT_air = self._dw_sat_dT(t_air_sat)
            d_omega_dtw = dw_dT_air * dT_air_dT_w
        else:
            denom = (h_sw - ha) + (le - 1.0) * ((h_sw - ha) - (w_sw - omega) * hv)
            d_omega_dtw = (cp_w * self.lg_ratio * (w_sw - omega)) / max(denom, 1e-7)
        return [d_ha_dtw, d_omega_dtw]

    def solve(self):
        """
        Integral from bottom to top (from t_out to t_in).
        Returns DataFrame with air profile along tower height.
        """
        ha_in = self.air_in.wet_air_enthalpy().magnitude
        omega_in = self.air_in.omega().magnitude
        y0 = [ha_in, omega_in]

        try:
            sol = solve_ivp(
                self.poppe_system,
                (self.t_out, self.t_in),
                y0,
                t_eval=np.linspace(self.t_out, self.t_in, 50),
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
                f"Reached tw={sol.t[-1]:.2f}°C of {self.t_in:.2f}°C"
            )
        air_temperatures, air_humidities = self.air_in.decode_results(sol.y[0], sol.y[1])
        zones = []
        fog_water = []
        for i, (h, w, t_a) in enumerate(zip(sol.y[0], sol.y[1], air_temperatures)):
            w_sw_air = psychrolib.GetHumRatioFromRelHum(t_a, 1.0, self.air_in.pressure)
            if w > w_sw_air:
                zones.append("fog")
                fog_water.append((w - w_sw_air))  
            else:
                zones.append("unsaturated")
                fog_water.append(0.0)

        df = pd.DataFrame({
            "water_temp_c":      sol.t,
            "air_temp_c":        air_temperatures,
            "air_rh_perc":       air_humidities,
            "air_omega_kg_kg":   sol.y[1],
            "air_enthalpy_j_kg": sol.y[0],
            "fog_kg_kg":          fog_water,   # туман в г воды на кг сухого воздуха
            "zone":              zones
        })
        omega_out = sol.y[1][-1]
        omega_in  = sol.y[1][0]

        # T_air на выходе
        t_air_out = air_temperatures[-1]
        w_sw_air_out = psychrolib.GetHumRatioFromRelHum(t_air_out, 1.0, self.air_in.pressure)

        if omega_out > w_sw_air_out:
            # Часть omega — это туман (капли), часть — пар
            omega_vapor_out = w_sw_air_out          # пар
            omega_fog_out   = omega_out - w_sw_air_out  # туман (капли)
        else:
            omega_vapor_out = omega_out
            omega_fog_out   = 0.0

        # Испарение и унос на 1 кг циркулирующей воды
        G_per_L = 1.0 / self.lg_ratio  # кг воздуха на кг воды

        evaporation  = (omega_vapor_out - omega_in) * G_per_L  # кг пара / кг воды
        fog_carryover = omega_fog_out * G_per_L                 # кг тумана / кг воды
        total_loss   = evaporation + fog_carryover              # кг воды / кг воды

        df.attrs['evaporation_kg_kg']   = evaporation
        df.attrs['fog_carryover_kg_kg'] = fog_carryover
        df.attrs['total_water_loss_kg_kg'] = total_loss
        h_sw_arr = [self.air_in.saturated_air_enthalpy(Q_(tw, u.degC)).magnitude 
                    for tw in df['water_temp_c']]
        cp_w_arr = [self.water_in.specific_heat(temperature=Q_(tw, u.degC)).magnitude 
                    for tw in df['water_temp_c']]
        integrand = [cp / max(h_sw - ha, 1.0) 
                    for cp, h_sw, ha in zip(cp_w_arr, h_sw_arr, df['air_enthalpy_j_kg'])]
        me = trapezoid(integrand, df['water_temp_c'])  
        df.attrs["merkel_number"] = me      
        return df
    @staticmethod
    def process_results(df):
        return dict(
            evaporation = df.attrs["total_water_loss_kg_kg"],
            me = df.attrs["merkel_number"],
            fog = df.attrs["fog_carryover_kg_kg"],
            fog_force = float((df["zone"] == "fog").sum()/df["zone"].size*100.0))
