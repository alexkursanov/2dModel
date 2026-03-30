from math import log, exp
import numpy as np
import scipy.optimize as optimize
from .constants import *


def calculate():
    def chi(v):
        if v <= 0.0:
            return chi_1 + chi_2 * v / v_max
        else:
            return chi_1

    def k_p_v(v):
        return chi(v) * chi_0 * q_v(v) * m_0 * G_star(v)

    def k_m_v(v):
        return chi_0 * q_v(v) * (1.0 - chi(v) * m_0 * G_star(v))

    def q_v(v):
        if v <= 0.0:
            return q_1 - q_2 * v / v_max
        elif (v <= v_st) and (0.0 < v):
            return (q_4 - q_3) * v / v_st + q_3
        else:
            return q_4 / (1.0 + beta_Q * (v - v_st) / v_max) ** (alpha_Q)

    #   v/v_max ch-to v
    def G_star(v):
        if v <= 0:
            return 1 + 0.6 * v / v_max
        elif v <= v_1:
            return P_star(v) / ((0.4 * a + 1.0) * v / (a * v_max) + 1.0)
        else:
            return (
                P_star(v)
                * exp(-alpha_G * ((v - v_1) / v_max) ** (alpha_P))
                / ((0.4 * a + 1.0) * v / (a * v_max) + 1.0)
            )

    #   v/v_max ch-to v
    def P_star(v):
        if v <= 0.0:
            return a * (1.0 + v / v_max) / (a - v / v_max)
        else:
            return (
                1.0
                + d_h
                - (d_h) ** (2.0)
                * a
                / (
                    a * d_h / gamma2 * (v / v_max) ** (2.0)
                    + (a + 1.0) * v / v_max
                    + a * d_h
                )
            )

    def n_1(l_1, version=2):
        if version == 0:
            if g_1 * l_1 + g_2 < 0.0:
                n_1 = 0.0
            elif g_1 * l_1 + g_2 < 1.0:
                n_1 = g_1 * l_1 + g_2
            else:
                n_1 = 1.0
        elif version == 2:
            if (g_1 * l_1 + g_2) * (
                n1_A + (n1_K - n1_A) / (n1_C + n1_Q * exp(-n1_B * l_1)) ** (1 / n1_nu)
            ) < 0.0:
                n_1 = 0.0
            elif (g_1 * l_1 + g_2) * (
                n1_A + (n1_K - n1_A) / (n1_C + n1_Q * exp(-n1_B * l_1)) ** (1 / n1_nu)
            ) < 1.0:
                n_1 = (g_1 * l_1 + g_2) * (
                    n1_A
                    + (n1_K - n1_A) / (n1_C + n1_Q * exp(-n1_B * l_1)) ** (1 / n1_nu)
                )
            else:
                n_1 = 1.0
        elif version == 1:
            r_1 = 2.25
            k_1 = 8.0

            if (g_1 * l_1 + g_2) * (1 - (r_1 * (0.5 - l_1)) ** (2 * k_1)) < 0.0:
                n_1 = 0.0
            elif (g_1 * l_1 + g_2) * (1 - (r_1 * (0.5 - l_1)) ** (2 * k_1)) < 1.0:
                n_1 = (g_1 * l_1 + g_2) * (1 - (r_1 * (0.5 - l_1)) ** (2 * k_1))
            else:
                n_1 = 1.0
        return n_1

    def M(A):
        return (
            (A / A_tot) ** (mu)
            * (1.0 + (k_mu) ** (mu))
            / ((A / A_tot) ** (mu) + (k_mu) ** (mu))
        )

    def L_oz(l_1):
        if l_1 <= s055:
            return (l_1 + S_0) / (s046 + S_0)
        else:
            return (S_0 + s055) / (s046 + S_0)

    def N0(l_0):
        return (r0 - beta_2 * (exp(alpha_2 * l_0) - 1)) / llambda

    def fi(l_1):
        return k_p_v(0) * M(A0) * n_1(l_1, 2) * L_oz(l_1) * (1.0 - N0(l_1)) - k_m_v(
            0
        ) * N0(l_1)

    def L(l_0):
        return (
            l_0
            + (log(beta_1) - log(r0 + beta_1 - beta_2 * (exp(alpha_2 * l_0) - 1)))
            / alpha_1
        )

    def pi_N_A(l_0, A, N):
        N_A = A_tot * s_c * N / (A * 1)
        # N_A = A_tot * s_c * N / (A * 1)
        if N_A < 0.0:
            pi_N_A = 1
        elif N_A <= 1.0:
            pi_N_A = (pi_min) ** (N_A)
        else:
            pi_N_A = pi_min
        return pi_N_A

    def delenie():
        l200 = log((r0 + beta_2) / beta_2) / alpha_2
        l100 = l200
        a = 0.9 * l100
        b = l100
        if fi(a) == 0.0:
            return a
        if fi(b) == 0.0:
            return b

        x = a + (b - a) / 2.0
        while abs(fi(x)) >= 0.0000001:
            x = a + (b - a) / 2.0
            if fi(x) < 0:
                a = x
            if fi(x) >= 0:
                b = x
        return x

    A0 = 0.00128
    A0 = 2.55e-6
    A0 = 0.00111979064656134
    print("A = ", A0)

    # part of code from N.Vikulova, function delenie
    l_2 = delenie()
    l_1 = L(l_2)
    N = N0(l_2)
    l_3 = log((r0 + beta_3) / beta_3) / alpha_3

    # l_1 = optimize.newton_krylov(fi, l_0)
    print("l_1 = ", l_1)
    # l_2 = L(l_1)
    print("l_2 = ", l_2)
    # N = N0(l_1)
    print("N = ", N)
    # l_3 = log((r0+beta_3)/beta_3)/alpha_3
    print("l_3 = ", l_3)

    # l_1 *=rest_length
    # l_2 *=rest_length
    # l_3 *=rest_length

    # N = 0.01
    A_off_1 = a_off * pi_N_A(l_1, A0, N) * exp(-k_A * A0)
    Ca_i = A_off_1 * A0 / (a_on * (A_tot - A0))
    print("Ca_i = ", Ca_i)

    l_0 = l_2 + l_3
    # l_0*=rest_length

    print("l_0 = ", l_0)
    return [A0, l_1, l_2, N, l_3, Ca_i], l_0


if __name__ == "__main__":
    calculate()
