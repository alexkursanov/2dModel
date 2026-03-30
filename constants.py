r0 = (
    2.55248904424517  # preload 2.55248904424517 - 100 %   0.733921 - 95% 0.199694 - 90%
)


g_CaL = 5.0e-5   # litre_per_farad_second( in L_type_Ca_current)1.
# *0.5
# *0.8
# *1.0
# *1.2
# *1.5
g_Ks = 0.392  # nanoS_per_picoF( in slow_time_dependent_potassium_current)CHANGED 0.8
# *0.5
# *0.8
# *1.0
# *1.2
# *1.5
g_Kr = 0.153  # nanoS_per_picoF( in rapid_time_dependent_potassium_current)CHANGED 0.8
# *0.5
# *0.8
# *1.0
# *1.2
# *1.5
g_K1 = 5.405  # nanoS_per_picoF( in inward_rectifier_potassium_current)
# *0.5
# *0.8
# *1.0
# *1.2
# *1.5

Vmax_up = 0.00765 # millimolar_per_millisecond( in calcium_dynamics)
# *0.8
# *1.0
# *1.2

g_bca = 0.000592  # nanoS_per_picoF( in calcium_background_current)
Buf_c = 0.13  # millimolar( in calcium_dynamics)
Buf_sr = 10.0  # millimolar( in calcium_dynamics)
Buf_ss = 0.4  # millimolar( in calcium_dynamics)
Ca_o = 2.0  # millimolar( in calcium_dynamics)
EC = 1.5  # millimolar( in calcium_dynamics)
K_buf_c = 0.00085  # millimolar( in calcium_dynamics)
K_buf_sr = 0.3  # millimolar( in calcium_dynamics)
K_buf_ss = 0.00025  # millimolar( in calcium_dynamics)
K_up = 0.00025  # millimolar( in calcium_dynamics)
V_leak = 0.00036  # per_millisecond( in calcium_dynamics)
V_rel = 0.1224  # per_millisecond( in calcium_dynamics)
V_sr = 0.001094  # micrometre3( in calcium_dynamics)
V_ss = 0.00005468  # micrometre3( in calcium_dynamics)
V_xfer = 0.00456  # per_millisecond( in calcium_dynamics)

k1_prime = 0.15  # per_millimolar2_per_millisecond( in calcium_dynamics)
k2_prime = 0.045  # per_millimolar_per_millisecond( in calcium_dynamics)
k3 = 0.06  # per_millisecond( in calcium_dynamics)
k4 = 0.005  # per_millisecond( in calcium_dynamics)
max_sr = 2.5  # dimensionless( in calcium_dynamics)
min_sr = 1.0  # dimensionless( in calcium_dynamics)
#
K_pCa = 0.0005  # millimolar( in calcium_pump_current)
g_pCa = 0.2476  # picoA_per_picoF( in calcium_pump_current)
g_Na = 14.838  # nanoS_per_picoF( in fast_sodium_current)CHANGED
Cm = 0.185  # microF( in membrane)
F = 96485.3415  # coulomb_per_millimole( in membrane)
R = 8314.472  # joule_per_mole_kelvin( in membrane)
T = 310.0  # kelvin( in membrane)
V_c = 0.016404  # micrometre3( in membrane)
stim_amplitude = 52.0  # picoA_per_picoF( in membrane)
stim_duration = 1.0  # millisecond( in membrane)
stim_period = 1000.0  # millisecond( in membrane)
stim_start = 10.0  # millisecond( in membrane)
K_o = 5.4  # millimolar( in potassium_dynamics)
g_pK = 0.0146  # nanoS_per_picoF( in potassium_pump_current)
g_Kr = 0.153  # nanoS_per_picoF( in rapid_time_dependent_potassium_current)CHANGED 0.8
P_kna = 0.03  # dimensionless( in reversal_potentials)

g_bna = 0.00029  # nanoS_per_picoF( in sodium_background_current)
K_NaCa = 10000.0  # picoA_per_picoF( in sodium_calcium_exchanger_current) CHANGED
K_sat = 0.1  # dimensionless( in sodium_calcium_exchanger_current)
Km_Ca = 1.38  # millimolar( in sodium_calcium_exchanger_current)
Km_Nai = 87.5  # millimolar( in sodium_calcium_exchanger_current)
# alpha = 2.5  # dimensionless( in sodium_calcium_exchanger_current)
alpha = 1.0  # dimensionless( in sodium_calcium_exchanger_current)CHANGED
gamma = 0.35  # dimensionless( in sodium_calcium_exchanger_current)
Na_o = 140.0  # millimolar( in sodium_dynamics)
K_mNa = 40.0  # millimolar( in sodium_potassium_pump_current)40
# K_mNa = 80.0  # millimolar( in sodium_potassium_pump_current) CHANGED
K_mk = 1.0  # millimolar( in sodium_potassium_pump_current)
P_NaK = 2.724  # picoA_per_picoF( in sodium_potassium_pump_current)
g_to = 0.294  # nanoS_per_picoF( in transient_outward_current) 2.5
stim_end = 1e15

Ca_sense = 0.0  # CHANGED from Vikulova TNNP
k_Ca_sense = 0.0069  # CHANGED from Vikulova TNNP
K_inh = 1.0

# ---------------------SR original set------------------------------
# Max_SR = 15
# Min_SR = 1
# EC_50_SR = 0.45
# H = 2.5
# k_o_Ca = 10
# k_i_Ca = 0.5
# k_im = 0.005
# k_om = 0.06
# k_s = 1.102  #25
##---------------------SR adaptiv for TP------------------------------
Max_SR = 2.5
Min_SR = 1.0
EC_50_SR = 1.5
H = 2.0
k_o_Ca = 2.1  # 0.15 #2.1
k_i_Ca = 0.025  # 0.045 #0.025
k_om = 0.06
k_im = 0.005
k_s = 0.1224
# ---------------------SR original set------------------------------
# Max_SR = 2.5
# Min_SR = 1.0
# EC_50_SR = 1.5
# H = 2.0
# k_o_Ca = 10
# k_i_Ca = 0.5
# k_im = 0.005
# k_om = 0.06
# k_s = 0.1224


# fibroblast

percent_fibrosis = 0.5
count_fibroblasts = 4
Cm_fi = 6.3  # 0.0063 # nanoS_per_picoF
g_gap = 3.0  # nano Simens
I_NaK_fib_max = 2.002  # 0.002002#pA/pF
K_mK_fib = 1.0  # mmol/L
K_mNa_fib = 11.0  # mmol/L
B_fib = -200  # mV
V_rev_fib = -150  # mV
E_k_fib = -87  # mV
g_Kv_fib = 0.25  # nS per pF
g_bNa_fib = 0.0095  # nS/pF
g_K1_fib = 0.4822  # #nS/pF
K_o_fib = 5.4  # millimolar( in potassium_dynamics)
Na_o_fib = 140.0  # millimolar( in sodium_dynamics)
V_c_fib = 0.016404
g_max = 0.0014400144001440014
mech_sensitive_chanels = True
alpha_fib = 14.6  # per_micrometre (in parameters)
beta_fib = 0.009  # millinewton (in parameters)


# EKB PARAMETERS

alpha_1 = 14.6  # per_micrometre (in parameters)
beta_1 = 4.2  # millinewton (in parameters)
alpha_2 = 14.6  # per_micrometre (in parameters)
beta_2 = 0.009  # millinewton (in parameters)
alpha_3 = 55.0  # per_micrometre (in parameters)
beta_3 = 0.11  # millinewton (in parameters)
llambda = 250.0  # millinewton (in parameters) 350
q_1 = 0.0173  # per_second (in parameters_izakov_et_al_1991)
q_2 = 0.259  # per_second (in parameters_izakov_et_al_1991)
q_3 = 0.0173  # per_second (in parameters_izakov_et_al_1991)
q_4 = 0.015  # per_second (in parameters_izakov_et_al_1991)
v_max = 0.0055  #  micrometre_per_second (in parameters)
a = 0.25  # dimensionless (in parameters)
alpha_Q = 10.0  # dimensionless (in parameters_izakov_et_al_1991)
beta_Q = 5.0  # dimensionless (in parameters_izakov_et_al_1991)
x_st = 0.964285  # dimensionless (in parameters_izakov_et_al_1991)
alpha_G = 1.0  # dimensionless (in parameters_izakov_et_al_1991)
m_0 = 0.9  # dimensionless (in parameters)
g_1 = 0.6  # per_micrometre (in crossbridge_kinetics)
g_2 = 0.52  # dimensionless (in crossbridge_kinetics)
S_0 = 1.14  # micrometre (in parameters_izakov_et_al_1991)
chi_1 = 0.55  # dimensionless (in parameters)
pi_min = 0.02
chi_2 = 0.0
d_h = 0.5  # dimensionless (in parameters)
m = 1.7  # WHATA
chi_0 = 2.1  # dimensionless (in parameters)
alpha_P = 4.0  # dimensionless (in parameters)
q_st = 1000.0  # WHATA
alpha_vp_l = 16.0  # per_micrometre (in CE_velocity)
alpha_vp_s = 16.0  # per_micrometre (in CE_velocity)
beta_vp_l = 0.1  # millinewton_second_per_micrometre (in CE_velocity)
beta_vp_s = 10  # millinewton_second_per_micrometre (in CE_velocity)
alpha_vs_l = 46.0  # per_micrometre (in PE_velocity) alp_vs
alpha_vs_s = 39.0  # per_micrometre (in PE_velocity) alp_vsr
beta_vs_l = 20.0  # millinewton_second_per_micrometre (in PE_velocity) vs2
beta_vs_s = 60.0  # millinewton_second_per_micrometre (in PE_velocity) vs2rel
a_off = 0.17  # per_second (in intracellular_calcium_concentration) #changed
a_on = 35.0
B_1_tot = 0.0  # millimolar (in intracellular_calcium_concentration)
B_2_tot = 0.0  # millimolar (in intracellular_calcium_concentration)
a_eqmin = 0.001299042
tau_inf = 1500
k_mu = 0.6  # dimensionless (in parameters)
mu = 3.3  # dimensionless (in parameters)
s_c = 1.0
k_A = 28.0  # per_millimolar (in parameters_izakov_et_al_1991) CHANGED
A_tot = 0.07  # millimolar (in intracellular_calcium_concentration)

rest_length = 1.67

n1_A = 0.5
n1_B = 55.0
n1_C = 1.0
n1_Q = 0.835
n1_K = 1.0
n1_nu = 5.0

s055 = 0.55
s046 = 0.46

v_st = x_st * v_max
v_1 = v_max / 10.0
gamma2 = a * d_h * (0.1) ** (2.0) / (3.0 * a * d_h - (a + 1.0) * 0.1)
F_afterload = 0.0

from initial_conditions import calculate as __calculate

l_0 = __calculate()[-1]
# l_0 = 0.2660143818868154*1.67
print("l_0 = ", l_0)
