cdef extern from "<math.h>" nogil:
    double pi "M_PI"  # as in Python's math module
    double exp(double x)
    double sqrt(double x)

import numpy as np
cimport numpy as np
from scipy.special.cython_special cimport erf
from scipy import optimize


cdef double phi_cdf(double x):
    #'Cumulative distribution function for the standard normal distribution'
    return (1.0 + erf(x / sqrt(2.0))) / 2.0


cdef double dgt(double g, double t):
    out = 1. - exp(-g * t)
    return out

cdef Pt_MKT_c(double time, np.ndarray[np.float64_t] rf_times,np.ndarray[np.float64_t] rf_values):
    cdef double  RateTime

    RateTime = np.interp(time, rf_times, rf_values)
    p_out = exp(- RateTime * time)

    return p_out


cdef price_swaption_c(np.ndarray[np.float64_t] prm_g2pp, double t_exp, double t_mat, double dt_d,double swp_atm_d,
                    np.ndarray[np.float64_t] rf_times,np.ndarray[np.float64_t] rf_values,int call_flag,int n_max):
    # definition
    cdef double t_mat_n, T
    cdef double a1, s1, b1, s2, rho1
    cdef double G2a, G2b, Gab, sigma_x, sigma_y, rho_xy, mu_x, mu_y

    cdef int k_max, k, i_max, i
    cdef double sum_int, x_min, x_max, dx, y_star_0, x_i
    cdef double y_star_n, h1_x, gauss_x, phi_h1_x, sum_phi, int_tmp, p_0_T, swpt_res

    t_mat_n = t_exp + t_mat
    T = t_exp

    a1 = prm_g2pp[0]
    s1 = prm_g2pp[1]
    b1 = prm_g2pp[2]
    s2 = prm_g2pp[3]
    rho1 = prm_g2pp[4]

    G2a = dgt(2. * a1, T)
    G2b = dgt(2. * b1, T)
    Gab = dgt(a1 + b1, T)

    sigma_x = s1 * np.sqrt(G2a / (2. * a1))
    sigma_y = s2 * np.sqrt(G2b / (2. * b1))
    rho_xy = (rho1 * s1 * s2) / ((a1 + b1) * sigma_x * sigma_y) * Gab

    mu_x = -M_x_T(0, T, T, a1, s1, b1, s2, rho1)
    mu_y = -M_y_T(0, T, T, a1, s1, b1, s2, rho1)

    # ---------- set coupon time -----------------

    k_max = int((t_mat_n - t_exp) / dt_d)

    cdef np.ndarray _coupon_time_array = np.zeros(k_max, dtype = 'double')
    cdef np.ndarray _coupon_array = np.zeros(k_max, dtype = 'double')

    for k in xrange(1, k_max):
        _coupon_time_array[k-1] = t_exp + dt_d * k
        _coupon_array[k-1] = swp_atm_d * dt_d

    _coupon_time_array[k_max-1] = t_exp + dt_d * k_max
    _coupon_array[k_max-1] = 1. + swp_atm_d * dt_d


    sum_int = 0.0

    x_min = mu_x - 5. * sigma_x
    x_max = mu_x + 5. * sigma_x

    dx = (x_max - x_min) / n_max

    i_max = int(n_max)

    y_star_0 = 0.0001
    for i in xrange(0, i_max):
        x_i = x_min + dx * i

        ff = optimize.root(compute_ystar, y_star_0, method='lm',
                args=(t_exp, x_i,  a1, s1, b1, s2, rho1,_coupon_time_array, _coupon_array, rf_times, rf_values))
        y_star_n = ff.x[0]

        #y_star_n = newtons_method(compute_ystar, compute_ystar_derivate, y_star_0,
        #                           t_exp, x_i, a1,s1,b1,s2,rho1, _coupon_time_array,
        #                           _coupon_array, rf_times, rf_values, 1e-12)

        h1_x = h1_function(x_i, y_star_n, mu_x, mu_y, sigma_x, sigma_y, rho_xy)

        gauss_x = exp(-0.5*(((x_i - mu_x)*(x_i - mu_x))/(sigma_x*sigma_x)))/(sigma_x*sqrt(2.*pi))

        phi_h1_x = phi_cdf(-call_flag * h1_x)

        sum_phi = compute_N(x_i, t_exp, t_mat_n, sigma_x, sigma_y, rho_xy, mu_x, mu_y, y_star_n, call_flag, a1, s1, b1, s2, rho1,
                            _coupon_array, _coupon_time_array, rf_times, rf_values)

        int_tmp = gauss_x * (phi_h1_x - sum_phi)

        sum_int = sum_int + int_tmp * dx

    p_0_T = price_model(0, t_exp, 0., 0., a1,s1,b1,s2,rho1, rf_times, rf_values)

    swpt_res = call_flag * p_0_T * sum_int

    return swpt_res


cdef double newtons_method(f, df,double x0, double a1,double a2,double a3, double a4, double a5, double a6, double a7,
                            np.ndarray[np.float64_t] a8,np.ndarray[np.float64_t] a9,np.ndarray[np.float64_t] a10,np.ndarray[np.float64_t] a11, double e):
    cdef double delta

    delta = abs(0. - f(x0, a1, a2, a3,a4,a5,a6,a7,a8,a9,a10,a11))
    while delta > e:
        x0 = x0 - f(x0, a1, a2, a3,a4,a5,a6,a7,a8,a9,a10,a11) / df(x0, a1, a2, a3,a4,a5,a6,a7,a8,a9,a10,a11)
        delta = abs(0. - f(x0, a1, a2, a3,a4,a5,a6,a7,a8,a9,a10,a11))

    return x0


cdef double h1_function(double x, double y_bar, double mu_x, double mu_y, double sigma_x, double sigma_y, double rho_xy):
    cdef double rde, z1, z2

    rde = np.sqrt(1. - rho_xy * rho_xy)
    z1 = (y_bar - mu_y) / (sigma_y * rde)
    z2 = rho_xy * (x - mu_x) / (sigma_x * rde)

    return z1 - z2


cdef compute_ystar(double y_star, double _t_exp,double _x, double a1,double s1,double b1,double s2,double rho1,
                            _coupon_time_array, _coupon_array, rf_times, rf_values):
    cdef double tmp, T, PtT, cpn_df_tmp, aa
    cdef int j, jmax

    tmp = 0.0
    T = _t_exp

    jmax = len(_coupon_time_array)

    for j in xrange(0, jmax):
        PtT = price_model(T, _coupon_time_array[j], _x, y_star, a1,s1,b1,s2,rho1, rf_times, rf_values)

        cpn_df_tmp = _coupon_array[j] * PtT
        tmp = tmp + cpn_df_tmp

    aa = tmp - 1.0

    return aa


cdef compute_ystar_derivate(double y_star,double _t_exp,double _x, double a1,double s1,double b1,double s2,double rho1,
                            _coupon_time_array, _coupon_array, rf_times, rf_values):
    cdef double tmp, T, PtT, cpn_df_tmp, aa, t_i
    cdef int j, jmax

    tmp = 0.0
    T = _t_exp

    jmax = len(_coupon_time_array)

    for j in xrange(0, jmax):
        t_i = _coupon_time_array[j]

        PtT = price_model(T, t_i, _x, y_star, a1,s1,b1,s2,rho1, rf_times, rf_values)

        cpn_df_tmp = _coupon_array[j] * PtT * (-B_ztT(b1, T, t_i))
        tmp = tmp + cpn_df_tmp

    aa = tmp

    return aa


cdef double M_x_T(double s,double t,double T,double a1,double s1,double b1,double s2,double rho1):
    cdef double g1a_Tt, g1a_Tt2s, g1b_Tt, g1ab, G1, G2, G3, h1, h2, h3, m_x_T

    g1a_Tt = exp(-a1 * (T - t))
    g1a_Tt2s = exp(-a1 * (T + t - 2. * s))

    g1b_Tt = exp(-b1 * (T - t))
    g1ab = exp(-b1 * T - a1 * t + (a1 + b1) * s)

    G1 = dgt(a1, t - s)
    G2 = g1a_Tt - g1a_Tt2s
    G3 = g1b_Tt - g1ab

    h1 = (s1 * s1) / (a1 * a1) + (rho1 * s1 * s2) / (a1 * b1)
    h2 = (s1 * s1) / (2.0 * a1 * a1)
    h3 = (rho1 * s1 * s2) / (b1 * (a1 + b1))

    m_x_T = h1 * G1 - h2 * G2 - h3 * G3

    return m_x_T


cdef double M_y_T(double s,double t,double T,double a1,double s1,double b1,double s2,double rho1):
    cdef double g1b_Tt, g1b_Tt2s, g1a_Tt, g1ba, G1, G2, G3, h1, h2, h3, m_y_T

    g1b_Tt = exp(-b1 * (T - t))
    g1b_Tt2s = exp(-b1 * (T + t - 2. * s))

    g1a_Tt = exp(-a1 * (T - t))
    g1ba = exp(-a1 * T - b1 * t + (a1 + b1) * s)

    G1 = dgt(b1, t - s)
    G2 = g1b_Tt - g1b_Tt2s
    G3 = g1a_Tt - g1ba

    h1 = (s2 * s2) / (b1 * b1) + (rho1 * s1 * s2) / (a1 * b1)
    h2 = (s2 * s2) / (2. * b1 * b1)
    h3 = (rho1 * s1 * s2) / (a1 * (a1 + b1))

    m_y_T = h1 * G1 - h2 * G2 - h3 * G3

    return m_y_T


cdef double compute_N(double x, double t_exp, double t_mat, double sigma_x, double sigma_y, double rho_xy, double mu_x, double mu_y, double y_star,
              int call_flag, double a1,double s1,double b1,double s2,double rho1,
              np.ndarray[np.float64_t] _coupon_array, np.ndarray[np.float64_t] _coupon_time_array, np.ndarray[np.float64_t] rf_times,np.ndarray[np.float64_t] rf_values):

    cdef double sum_n, T, c_i, t_i, b_i, lambda_x, h1_x, h2_x, k_x, phi_x
    cdef int j_max, j

    sum_n = 0.0

    T = t_exp

    j_max = len(_coupon_time_array)

    for j in xrange(0, j_max):
        c_i = _coupon_array[j]
        t_i = _coupon_time_array[j]

        b_i = B_ztT(b1, T, t_i)

        lambda_x = c_i * A2_tT(T, t_i, rf_times, rf_values, a1, s1, b1, s2, rho1) * exp(-B_ztT(a1, T, t_i) * x)

        h1_x = h1_function(x, y_star, mu_x, mu_y, sigma_x, sigma_y, rho_xy)

        h2_x = h1_x + b_i * sigma_y * np.sqrt(1. - rho_xy * rho_xy)

        k_x = - b_i * (mu_y - 0.5 * (1. - rho_xy * rho_xy) * sigma_y * sigma_y * b_i + rho_xy * sigma_y * (
                (x - mu_x) / sigma_x))

        phi_x = phi_cdf(-call_flag * h2_x)

        sum_n = sum_n + lambda_x * exp(k_x) * phi_x

    return sum_n


cdef double B_ztT(double z,double t,double T):
    cdef double BtT
    BtT = dgt(z, T - t) / z

    return BtT


cdef double A2_tT(double t, double T, np.ndarray[np.float64_t] rf_times, np.ndarray[np.float64_t] rf_values,
                    double a1,double s1,double b1,double s2,double rho1):
    cdef double P_mkt_0_t,P_mkt_0_T,v_tT,v_0T,v_0t,A1,AtT

    P_mkt_0_t = Pt_MKT_c(t, rf_times, rf_values)
    P_mkt_0_T = Pt_MKT_c(T, rf_times, rf_values)

    v_tT = V_tT(t, T, a1,s1,b1,s2,rho1)
    v_0T = V_tT(0, T, a1,s1,b1,s2,rho1)
    v_0t = V_tT(0, t, a1,s1,b1,s2,rho1)

    A1 = 0.5 * (v_tT - v_0T + v_0t)
    AtT = (P_mkt_0_T / P_mkt_0_t) * exp(A1)

    return AtT


cdef double V_tT(double t, double T, double a,double sigma,double b,double eta,double rho):
    cdef double z1, z2, z3

    z1 = (sigma * sigma) / (a * a) * (
            T - t + (2 / a) * exp(-a * (T - t)) - (1 / (2 * a)) * exp(-2 * a * (T - t)) - (3 / (2 * a)))
    z2 = (eta * eta) / (b * b) * (
            T - t + (2 / b) * exp(-b * (T - t)) - (1 / (2 * b)) * exp(-2 * b * (T - t)) - (3 / (2 * b)))
    z3 = 2 * rho * (sigma * eta) / (a * b) * (T - t + (exp(-a * (T - t)) - 1) / a +
                                              (exp(-b * (T - t)) - 1) / b - (exp(-(a + b) * (T - t)) - 1) / (
                                                          a + b))
    return z1 + z2 + z3


cdef double price_model(double t,double T,double xt, double yt, double a1,double s1,double b1,double s2,double rho1,np.ndarray[np.float64_t] rf_times,np.ndarray[np.float64_t] rf_values):
    cdef double a, Ba,Bb, out

    a = A2_tT(t, T, rf_times, rf_values, a1, s1, b1, s2, rho1)
    Ba = B_ztT(a1, t, T)
    Bb = B_ztT(b1, t, T)
    out = a * exp(- Ba * xt - Bb * yt)

    return out


#### cdef double loss_g2pp_model(np.ndarray[np.float64_t] list_model_params,
####                         np.ndarray[np.float64_t] t_exp_vector,
####                         np.ndarray[np.float64_t] t_mat_vector,
####                         np.ndarray[np.float64_t] swp_atm_d_vector,
####                         np.ndarray[np.float64_t] mkt_price_vector,
####                         np.ndarray[np.float64_t] rf_times,
####                         np.ndarray[np.float64_t] rf_values,
####                         double tenor,
####                         int call_flag):

cdef double loss_g2pp_model(list_model_params,
                        t_exp_vector,
                        t_mat_vector,
                        swp_atm_d_vector,
                        mkt_price_vector,
                        rf_times,
                        rf_values,
                        double tenor,
                        int call_flag):

    cdef double diff_sum,t_exp,t_mat,swp_atm_d,model_price_tmp,mkt_price_tmp,diff
    cdef int i, i_max

    diff_sum = 0.0
    i_max = int(len(t_exp_vector))

    for i in xrange(0,i_max):
        t_exp = t_exp_vector[i]
        t_mat = t_mat_vector[i]
        swp_atm_d = swp_atm_d_vector[i]
        model_price_tmp = price_swaption_c(list_model_params, t_exp, t_mat, tenor, swp_atm_d, rf_times, rf_values, call_flag, n_max = 10)
        mkt_price_tmp = mkt_price_vector[i]
        diff = abs(model_price_tmp - mkt_price_tmp) #/ mkt_price_tmp
        diff = diff * diff
        diff_sum = diff_sum + diff

    return diff_sum

cpdef found_opt( np.ndarray[np.float64_t] InitialGuess,
               InitialBound,
               np.ndarray[np.float64_t] t_exp_vector,
               np.ndarray[np.float64_t] t_mat_vector,
               np.ndarray[np.float64_t] swp_atm_d_vector,
               np.ndarray[np.float64_t] mkt_price_vector,
               np.ndarray[np.float64_t] rf_times,
               np.ndarray[np.float64_t] rf_values,
               double tenor,
               int call_flag):

    ff = optimize.minimize(fun=loss_g2pp_model, method='TNC',
                            args=(t_exp_vector,t_mat_vector,swp_atm_d_vector,mkt_price_vector,rf_times,rf_values,tenor,call_flag),
                            x0= InitialGuess, bounds = InitialBound,
                            options={'disp': True})

    return ff

cpdef price_swaption(np.ndarray[np.float64_t] prm_g2pp, double t_exp, double t_mat, double dt_d,double swp_atm_d,
                        np.ndarray[np.float64_t] rf_times,np.ndarray[np.float64_t] rf_values,int call_flag,int n_max):
    return price_swaption_c(prm_g2pp, t_exp, t_mat, dt_d, swp_atm_d, rf_times, rf_values, call_flag, n_max)


cpdef Pt_MKT(double time, np.ndarray[np.float64_t] rf_times,np.ndarray[np.float64_t] rf_values):
    return Pt_MKT_c(time, rf_times, rf_values)

