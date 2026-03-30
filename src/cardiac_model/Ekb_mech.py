
import numpy as np
from numba import jit
from numba.experimental import jitclass
from .constants import *
from numba import float64, float32, boolean
from math import exp, log, floor, sqrt

spec = [
    ('t0', float32),
    ('F_afterload', float64),
    ('currents', float64[:]),
    ('new_coop', boolean)
]


@jitclass(spec)
class ETNNP(object):
    def __init__(self):
        self.t0 = 0.0
        self.F_afterload = 0.0
        self.currents = np.zeros(11, dtype=np.float64)
        self.new_coop = False

    @property
    def show_F(self):
        return self.F_afterload

    @property
    def return_currents(self):
        return self.currents

    def SetNewCoopOn(self):
        self.new_coop = True
    
    def SetNewCoopOff(self):
        self.new_coop = False 

    # EKATERINBURG FUNCTIONS

    def chi(self, v):
        if v <= 0.0:
            return chi_1 + chi_2 * v / v_max
        else:
            return chi_1

    def k_p_v(self, v):
        return self.chi(v) * chi_0 * self.q_v(v) * m_0 * self.G_star(v)

    def k_m_v(self, v):
        return chi_0 * self.q_v(v) * (1.0 - self.chi(v) * m_0 * self.G_star(v))

    def q_v(self, v):
        if (v <= 0.0):
            return q_1 - q_2 * v / v_max
        elif ((v <= v_st) and (0.0 < v)):
            return (q_4 - q_3) * v / v_st + q_3
        else:
            return q_4 / (1.0 + beta_Q * (v - v_st) / v_max)**alpha_Q

    def G_star(self, v):
        if v <= 0:
            return 1+0.6*v/v_max
        elif v <= v_1:
            return self.P_star(v) / ((0.4 * a + 1.0) * v / (a * v_max) + 1.0)
        else:
            return self.P_star(v) * exp(-alpha_G * ((v - v_1) / v_max)**alpha_P) / ((0.4 * a + 1.0) * v / (a * v_max) + 1.0)

    def P_star(self, v):
        if (v <= 0.0):
            return a * (1.0 + v / v_max) / (a - v / v_max)
        else:
            return 1.0 + d_h - (d_h**2.0) * a / (a * d_h / gamma2 * (v / v_max)**2.0 + (a + 1.0) * v / v_max + a * d_h)

    def n_1(self, l_1, version=2):
        if version == 0:
            if (g_1 * l_1 + g_2 < 0.0):
                n_1 = 0.0
            elif (g_1 * l_1 + g_2 < 1.0):
                n_1 = g_1 * l_1 + g_2
            else:
                n_1 = 1.0
        elif version == 2:
            if (g_1 * l_1 + g_2) * (n1_A + (n1_K - n1_A) / (n1_C + n1_Q * exp(-n1_B * l_1))**(1 / n1_nu)) < 0.0:
                n_1 = 0.0
            elif (g_1 * l_1 + g_2) * (n1_A + (n1_K - n1_A) / (n1_C + n1_Q * exp(-n1_B * l_1))**(1 / n1_nu)) < 1.0:
                n_1 = (g_1 * l_1 + g_2) * (n1_A + (n1_K - n1_A) / (n1_C + n1_Q * exp(-n1_B * l_1))**(1 / n1_nu))
            else:
                n_1 = 1.
        elif version == 1:
            r_1 = 2.25/rest_length
            k_1 = 8.

            if ((g_1 * l_1 + g_2) * (1 - (r_1 * (0.5 - l_1)) ** (2 * k_1)) < 0.0):
                n_1 = 0.0
            elif (g_1 * l_1 + g_2) * (1 - (r_1 * (0.5 - l_1)) ** (2 * k_1)) < 1.0:
                n_1 = (g_1 * l_1 + g_2) * (1 - (r_1 * (0.5 - l_1)) ** (2 * k_1))
            else:
                n_1 = 1.0
        return n_1

    def M(self, A):
        return (A / A_tot)**mu * (1.0 + (k_mu**mu)) / ((A / A_tot)**mu + (k_mu**mu))

    def L_oz(self, l_1):
        if (l_1 > s055):
            return (l_1 + S_0) / (s046 + S_0)
        else:
            return (S_0 + s055) / (s046 + S_0)
        # return 1

    def pi_N_A(self, l_1, N, A):
        N_A = A_tot * s_c * N / (A * self.L_oz(l_1))
        # N_A = A_tot * s_c * N / (A * 1)
    
        if (N_A <= 0.0):
            pi_N_A = 1.0
        elif (N_A <= 1.0):
            pi_N_A = (pi_min**N_A)
        else:
            pi_N_A = pi_min
        return pi_N_A


    def pi_N(self,l_1, N, A):
        if N < 0:
            pi = 1.
        elif A_tot*N >=0 and A_tot*N <= 1/s_c:
            pi = pi_min**(A_tot*N*s_c)
        elif A_tot*N >= 1/s_c:
            pi = pi_min
        return pi

    def main(self, time, Y):
        # TNNP variables
        d = Y[0]
        f2 = Y[1]
        fCass = Y[2]
        f = Y[3]
        Ca_SR = Y[4]
        Ca_i = Y[5]
        Ca_ss = Y[6]
        R_prime = Y[7]
        h = Y[8]
        j = Y[9]
        m = Y[10]
        V = Y[11]
        K_i = Y[12]
        Xr1 = Y[13]
        Xr2 = Y[14]
        Xs = Y[15]
        Na_i = Y[16]
        r = Y[17]
        s = Y[18]

        # Ekb variables

        v = Y[19]
        w = Y[20]
        N = Y[21]
        A = Y[22]
        l_1 = Y[23]
        l_2 = Y[24]
        l_3 = Y[25]

        machine_zero = 1e-15

        if self.F_afterload <= machine_zero:
            isotonic = False
        else:
            isotonic = True
        physiol_mode = False  # dimensionless ( in physiological)

        # ------------------------------------------------------------------------------
        # Computation
        # ------------------------------------------------------------------------------

        # Ekb
        if (v <= 0.0):
            alpha_p = alpha_vp_l
        else:
            alpha_p = alpha_vp_s

        if (v <= 0.0):
            k_P_vis = beta_vp_l * exp(alpha_vp_l * l_1)
        else:
            k_P_vis = beta_vp_s * exp(alpha_vp_s * l_1)

        K_chi = self.k_p_v(v) * self.M(A) * self.n_1(l_1, 2) * \
            self.L_oz(l_1) * (1.0 - N) - self.k_m_v(v) * N
        
        p_v = self.P_star(v) / self.G_star(v)

        case_1 = a * (0.4 + 0.4 * a) / (v_max * ((a + 1.0) * 0.4)**2.0)
        case_2 = a * 1.0 * (1.0 + 0.4 * a + 1.2 * v / v_max + 0.6 * (v / v_max)**2.0) / (
            v_max * ((a - v / v_max) * (1.0 + 0.6 * v / v_max))**2.0)
        case_3 = (0.4 * a + 1.0) / (a * v_max)
        case_4 = 1.0 / v_max * exp(alpha_G * (v / v_max - v_1 / v_max)**alpha_P) * (
            (0.4 * a + 1.0) / a + alpha_G * alpha_P * (1.0 + (0.4 * a + 1.0) * v / (a * v_max)) * (v / v_max - v_1 / v_max)**(alpha_P - 1.0))

        if (v <= -v_max):
            p_prime_v = case_1
        elif ((-v_max < v) and (v <= 0.0)):
            p_prime_v = case_2
        elif ((0.0 < v) and (v <= v_1)):
            p_prime_v = case_3
        else:
            p_prime_v = case_4

        F_XSE = beta_3 * (exp(alpha_3 * l_3) - 1.0)
        F_muscle = F_XSE
        l = l_2 + l_3

        #CHANGED
        if ((isotonic) and (F_muscle > self.F_afterload) and (l <= l_0 * (1.0 + 1.0e-4))): 
            isotonic_mode = True
        else:
            isotonic_mode = False

        if ((time - floor(time / stim_period) * stim_period >= stim_start) and  (time - floor(time / stim_period) * stim_period <= stim_start + stim_duration)):
            i_Stim = -stim_amplitude
        else:
            i_Stim = 0.0

        if (isotonic_mode):
            phi_chi = -(llambda * K_chi * p_v + alpha_p * k_P_vis * v**2.0 + alpha_2 * beta_2 * exp(alpha_2 * l_2) * w) / (llambda * N * p_prime_v + k_P_vis)
        else:
            phi_chi = -(llambda * K_chi * p_v + alpha_p * k_P_vis * v**2.0 + (alpha_2 * beta_2 * exp(alpha_2 * l_2) + alpha_3 * beta_3 * exp(alpha_3 * l_3)) * w) / (llambda * N * p_prime_v + k_P_vis)

        if (isotonic_mode):
            phi_chi_2 = alpha_1 * beta_1 * exp(alpha_1 * (l_2 - l_1)) * v / (alpha_1 * beta_1 * exp(alpha_1 * (l_2 - l_1)) + alpha_2 * beta_2 * exp(alpha_2 * l_2))
        else:
            phi_chi_2 = alpha_1 * beta_1 * exp(alpha_1 * (l_2 - l_1)) * v / (alpha_1 * beta_1 * exp(alpha_1 * (l_2 - l_1)) + alpha_2 * beta_2 * exp(alpha_2 * l_2) + alpha_3 * beta_3 * exp(alpha_3 * l_3))

        if (w <= v):
            k_S_vis = beta_vs_l * exp(alpha_vs_l * (l_2 - l_1))
        else:
            k_S_vis = beta_vs_s * exp(alpha_vs_s * (l_2 - l_1))

        if (-machine_zero <= k_S_vis <= machine_zero):
            dv = (alpha_1 * beta_1 * exp(alpha_1 * (l_2 - l_1)) * (phi_chi_2 - v) - (llambda * K_chi * p_v + alpha_p * k_P_vis * v**2.0)) / (llambda * N * p_prime_v + k_P_vis)
        else:
            dv = phi_chi


        if (w <= v):
            alp_s = alpha_vs_l
        else:
            alp_s = alpha_vs_s

        if (isotonic_mode and abs(k_S_vis) > machine_zero):
            dw = (k_S_vis * (phi_chi - alp_s * (w - v)**(2.0)) - alpha_1 * beta_1 * exp(alpha_1 * (l_2 - l_1)) * (w - v) - alpha_2 * beta_2 * exp(alpha_2 * l_2) * w) / k_S_vis
        elif (not isotonic_mode and abs(k_S_vis) > machine_zero):
            dw = phi_chi - alp_s * (w - v)**2.0 - (alpha_1 * beta_1 * exp(alpha_1 * (l_2 - l_1)) * (w - v) + (alpha_2 * beta_2 * exp(alpha_2 * l_2) + alpha_3 * beta_3 * exp(alpha_3 * l_3)) * w) / k_S_vis
        elif abs(k_S_vis) <= machine_zero:
            dw = 0.0

        dN = K_chi

        F_CE = llambda * p_v * N
        F_SE = beta_1 * (exp(alpha_1 * (l_2 - l_1)) - 1.0)
        F_PE = beta_2 * (exp(alpha_2 * l_2) - 1.0)
        # F = F_SE + F_PE


        A_off_1 = a_off * self.pi_N_A(l_1, N, A) * exp(-k_A * A)
        A_off_2 = a_on * a_eqmin

        if self.new_coop:
                
            # Alpha = 1.0
            if A_off_1 > A_off_2:
                self.t0 = time
                Alpha = 1.0
            else:
                Alpha = exp((self.t0 - time) / tau_inf)
            A_off = Alpha * A_off_1 + (1.0 - Alpha) * A_off_2

            dA = a_on * (A_tot - A) * Ca_i - A_off * A
        else:
            dA = a_on * (A_tot - A) * Ca_i - A_off_1 * A



        dl_1_dt = v

        dl_1 = dl_1_dt

        if abs(k_S_vis) <= machine_zero:
            dl_2_dt = phi_chi_2
        else:
            dl_2_dt = w

        # if isotonic_mode:
        #     dl_2_dt = 0

        dl_2 = dl_2_dt

        
        if isotonic_mode:
            dl_3_dt = 0.0
        elif not isotonic_mode and abs(k_S_vis) <= machine_zero:
            dl_3_dt = -phi_chi_2
        elif not isotonic_mode and abs(k_S_vis) > machine_zero:
            dl_3_dt = -w

        if ((time >= stim_start) and (time <= stim_end) and (time - stim_start - floor((time - stim_start) / stim_period) * stim_period <= stim_duration)):
            dl_3 = -(l_2 + l_3 - l_0) / stim_duration
        else:
            dl_3 = dl_3_dt
        # dl_3 = dl_3_dt
        # TNNP

        i_CaL = g_CaL * d * f * f2 * fCass * 4.0 * (V - 15.0) * F**2.0 / (R * T) * (
            0.25 * Ca_ss * exp(2.0 * (V - 15.0) * F / (R * T)) - Ca_o) / (exp(2.0 * (V - 15.0) * F / (R * T)) - 1.0)
        d_inf = 1.0 / (1.0 + exp((-8.0 - V) / 7.5))
        alpha_d = 1.4 / (1.0 + exp((-35.0 - V) / 13.0)) + 0.25
        beta_d = 1.4 / (1.0 + exp((V + 5.0) / 5.0))
        gamma_d = 1.0 / (1.0 + exp((50.0 - V) / 20.0))
        tau_d = 1.0 * alpha_d * beta_d + gamma_d
        dd = (d_inf - d) / tau_d
        f2_inf = 0.67 / (1.0 + exp((V + 35.0) / 7.0)) + 0.33
        tau_f2 = 562.0 * exp(-(V + 27.0)**2.0 / 240.0) + 31.0 / (1.0 + exp((25.0 - V) / 10.0)) + 80.0 / (
            1.0 + exp((V + 30.0) / 10.0))
        df2 = (f2_inf - f2) / tau_f2
        fCass_inf = 0.6 / (1.0 + (Ca_ss / 0.05)**2.0) + 0.4
        tau_fCass = 80.0 / (1.0 + (Ca_ss / 0.05)**2.0) + 2.0
        dfCass = (fCass_inf - fCass) / tau_fCass
        f_inf = 1.0 / (1.0 + exp((V + 20.0) / 7.0))
        tau_f = 1102.5 * exp(-(V + 27.0)**(2.0) / 225.0) + 200.0 / (1.0 + exp((13.0 - V) / 10.0)) + 180.0 / (
            1.0 + exp((V + 30.0) / 10.0)) + 20.0
        df = (f_inf - f) / tau_f
        E_Ca = 0.5 * R * T / F * log(Ca_o / Ca_i)
        i_b_Ca = g_bca * (V - E_Ca)
        kcasr = max_sr - (max_sr - min_sr) / (1.0 + (EC / Ca_SR)**2.0)
        k1 = k1_prime / kcasr
        O = k1 * Ca_ss**2.0 * R_prime / (k3 + k1 * Ca_ss**2.0)
        i_rel = V_rel * O * (Ca_SR - Ca_ss)
        # CHANGED there was no inhibition and no K_inh 
        i_up = Vmax_up / (1.0 + K_up**2.0 / Ca_i**2.0)
        i_leak = V_leak * (Ca_SR - Ca_i)
        i_xfer = V_xfer * (Ca_ss - Ca_i)
        k2 = k2_prime * kcasr
        dR_prime = -k2 * Ca_ss * R_prime + k4 * (1.0 - R_prime)
        Ca_i_bufc = 1.0 / (1.0 + Buf_c * K_buf_c / (Ca_i + K_buf_c)**2.0)
        Ca_sr_bufsr = 1.0 / (1.0 + Buf_sr * K_buf_sr / (Ca_SR + K_buf_sr)**2.0)
        Ca_ss_bufss = 1.0 / (1.0 + Buf_ss * K_buf_ss / (Ca_ss + K_buf_ss)**2.0)
        i_p_Ca = g_pCa * Ca_i / (Ca_i + K_pCa)
        i_NaCa_first = exp(gamma * V * F / (R * T)) * (Na_i**3.0) * Ca_o
        i_NaCa_second = exp((gamma - 1.0) * V * F / (R * T)) * (Na_o**3.0) * Ca_i * alpha
        Ca_sense_factor = 1/(1+Ca_sense*Ca_i/k_Ca_sense) #CHANGED
        Sat_factor = 1/(1+K_sat*exp((gamma-1)*V*F/(R*T)))   # CHANGED
        #next is CHANGED, there are no Ca_sence_factor 
        i_NaCa = K_NaCa *(i_NaCa_first - i_NaCa_second) / (((Km_Nai**3.0) + (Na_o**3.0)) * (Km_Ca + Ca_o) * (1.0 + K_sat * exp((gamma - 1.0) * V * F / (R * T))))
        dCa_i = Ca_i_bufc * ((i_leak - i_up) * V_sr / V_c + i_xfer - 1.0 *(i_b_Ca + i_p_Ca - 2.0 * i_NaCa) * Cm / (2.0 * 1.0 * V_c * F) - dA)
        dCa_SR = Ca_sr_bufsr * (i_up - (i_rel + i_leak))
        dCa_ss = Ca_ss_bufss * (-1.0 * i_CaL * Cm / (2.0 * 1.0 * V_ss * F) + i_rel * V_sr / V_ss - i_xfer * V_c / V_ss)
        E_Na = R * T / F * log(Na_o / Na_i)
        i_Na = g_Na * (m**3.0) * h * j * (V - E_Na)
        h_inf = 1.0 / (1.0 + exp((V + 71.55) / 7.43))**2.0

        if (V < -40.0):
            alpha_h = 0.057 * exp(-(V + 80.0) / 6.8)
        else:
            alpha_h = 0.0

        if (V < -40.0):
            beta_h = 2.7 * exp(0.079 * V) + 310000.0 * exp(0.3485 * V)
        else:
            beta_h = 0.77 / (0.13 * (1.0 + exp((V + 10.66) / -11.1)))

        tau_h = 1.0 / (alpha_h + beta_h)
        dh = (h_inf - h) / tau_h
        j_inf = 1.0 / (1.0 + exp((V + 71.55) / 7.43))**2.0

        if (V < -40.0):
            alpha_j = (-25428.0 * exp(0.2444 * V) - 6.948e-6 * exp(-0.04391 * V)) * (V + 37.78) / 1.0 / (
                1.0 + exp(0.311 * (V + 79.23)))
        else:
            alpha_j = 0.0

        if (V < -40.0):
            beta_j = 0.02424 * exp(-0.01052 * V) / \
                (1.0 + exp(-0.1378 * (V + 40.14)))
        else:
            beta_j = 0.6 * exp(0.057 * V) / \
                (1.0 + exp(-0.1 * (V + 32.0)))

        tau_j = 1.0 / (alpha_j + beta_j)
        dj = (j_inf - j) / tau_j
        m_inf = 1.0 / (1.0 + exp((-56.86 - V) / 9.03))**2.0
        alpha_m = 1.0 / (1.0 + exp((-60.0 - V) / 5.0))
        beta_m = 0.1 / (1.0 + exp((V + 35.0) / 5.0)) + \
            0.1 / (1.0 + exp((V - 50.0) / 200.0))
        tau_m = 1.0 * alpha_m * beta_m
        dm = (m_inf - m) / tau_m
        E_K = R * T / F * log(K_o / K_i)
        alpha_K1 = 0.1 / (1.0 + exp(0.06 * (V - E_K - 200.0)))
        beta_K1 = (3.0 * exp(0.0002 * (V - E_K + 100.0)) + exp(0.1 * (V - E_K - 10.0))) / (
            1.0 + exp(-0.5 * (V - E_K)))
        xK1_inf = alpha_K1 / (alpha_K1 + beta_K1)
        i_K1 = g_K1 * xK1_inf * sqrt(K_o / 5.4) * (V - E_K)



        i_to = g_to * r * s * (V - E_K)
        i_Kr = g_Kr * sqrt(K_o / 5.4) * Xr1 * Xr2 * (V - E_K)
        E_Ks = R * T / F * log((K_o + P_kna * Na_o) / (K_i + P_kna * Na_i))
        i_Ks = g_Ks * (Xs**2.0) * (V - E_Ks)
        i_NaK = P_NaK * K_o / (K_o + K_mk) * Na_i / (Na_i + K_mNa) / (
            1.0 + 0.1245 * exp(-0.1 * V * F / (R * T)) + 0.0353 * exp(-V * F / (R * T)))
        i_b_Na = g_bna * (V - E_Na)
        i_p_K = g_pK * (V - E_K) / (1.0 + exp((25.0 - V) / 5.98))
        dV = -1.0 / 1.0 * (
            i_K1 + i_to + i_Kr + i_Ks + i_CaL + i_NaK + i_Na + i_b_Na + i_NaCa + i_b_Ca + i_p_K + i_p_Ca + i_Stim)
        dK_i = -1.0 * (i_K1 + i_to + i_Kr + i_Ks + i_p_K +
                       i_Stim - 2.0 * i_NaK) / (1.0 * V_c * F) * Cm
        xr1_inf = 1.0 / (1.0 + exp((-26.0 - V) / 7.0))
        alpha_xr1 = 450.0 / (1.0 + exp((-45.0 - V) / 10.0))
        beta_xr1 = 6.0 / (1.0 + exp((V + 30.0) / 11.5))
        tau_xr1 = 1.0 * alpha_xr1 * beta_xr1
        dXr1 = (xr1_inf - Xr1) / tau_xr1
        xr2_inf = 1.0 / (1.0 + exp((V + 88.0) / 24.0))
        alpha_xr2 = 3.0 / (1.0 + exp((-60.0 - V) / 20.0))
        beta_xr2 = 1.12 / (1.0 + exp((V - 60.0) / 20.0))
        tau_xr2 = 1.0 * alpha_xr2 * beta_xr2
        dXr2 = (xr2_inf - Xr2) / tau_xr2
        xs_inf = 1.0 / (1.0 + exp((-5.0 - V) / 14.0))
        alpha_xs = 1400.0 / sqrt(1.0 + exp((5.0 - V) / 6.0))
        beta_xs = 1.0 / (1.0 + exp((V - 35.0) / 15.0))
        tau_xs = 1.0 * alpha_xs * beta_xs + 80.0
        dXs = (xs_inf - Xs) / tau_xs
        dNa_i = -1.0 * (i_Na + i_b_Na + 3.0 * i_NaK +
                        3.0 * i_NaCa) / (1.0 * V_c * F) * Cm
        r_inf = 1.0 / (1.0 + exp((20.0 - V) / 6.0))
        tau_r = 9.5 * exp(-(V + 40.0)**2.0 / 1800.0) + 0.8
        dr = (r_inf - r) / tau_r
        s_inf = 1.0 / (1.0 + exp((V + 20.0) / 5.0))
        tau_s = 85.0 * exp(-(V + 45.0)**2.0 / 320.0) + \
            5.0 / (1.0 + exp((V - 20.0) / 5.0)) + 3.0
        ds = (s_inf - s) / tau_s
        # self.currents = np.array([I_Na, INaL,Ito,ICaL,ICaNa,ICaK,IKr,IKs,IK1, INaCa_i, INaCa_ss, INaK, INab, IKb, IpCa, ICab, Istim], dtype = np.float32)

        self.currents = np.array([i_Na, i_CaL, i_NaCa, i_K1, i_Kr, i_Ks, i_b_Ca, i_p_Ca, Ca_i_bufc, Ca_sr_bufsr, Ca_ss_bufss, k_P_vis, k_S_vis, p_v], dtype = np.float64)
        return np.array([dd, df2, dfCass, df, dCa_SR, dCa_i, dCa_ss, dR_prime, dh, dj, dm, dV, dK_i, dXr1, dXr2, dXs, dNa_i, dr, ds, dv, dw, dN, dA, dl_1, dl_2, dl_3])


if __name__ == "__main__":
    obj = ETNNP()

    Y0 = np.zeros([26])
    Y0[0] = 3.373e-5
    Y0[1] = 0.9755
    Y0[2] = 0.9953
    Y0[3] = 0.7888
    Y0[4] = 3.64
    Y0[5] = 0.000126
    Y0[6] = 0.00036
    Y0[7] = 0.9073
    Y0[8] = 0.7444
    Y0[9] = 0.7045
    Y0[10] = 0.00172
    Y0[11] = -85.23
    Y0[12] = 136.89
    Y0[13] = 0.00621
    Y0[14] = 0.4712
    Y0[15] = 0.0095
    Y0[16] = 8.604
    Y0[17] = 2.42e-8
    Y0[18] = 0.999998

    Y0[19] = 0.0  # v        (micrometre_per_second)( in CE_velocity)
    Y0[20] = 0.0  # w        (micrometre_per_second)( in PE_velocity)
    Y0[21] = 2.726318970e-6  # N
    Y0[22] = 6.7e-5  # A
    Y0[23] = 0.436321675  # l_1(micrometre)( in length)
    Y0[24] = 0.436328344  # l_2(micrometre)( in length)
    Y0[25] = 0.088805830  # l_3(micrometre)( in length)

    obj.main(0.0,  Y0)
