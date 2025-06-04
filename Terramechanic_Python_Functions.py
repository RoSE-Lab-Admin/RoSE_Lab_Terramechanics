import numpy as np
from scipy.integrate import quad

def calculate_r_s(r, h): #Effective Shearing Radius
  return r+0.5*h

def calculate_i(v, r_s, omega): #Wheel Slip
    if v < r_s * omega:
        return 1 - (v / (r_s * omega))
    elif v > r_s * omega:
        return (r_s * omega / v) - 1
    else:
        return 0

def calculate_theta_m(a1, a2, i, theta_1): #Angle of Maximum Stress
    return (a1 + a2 * i) * theta_1

def calculate_n(n0, n1, i): #Sinkage Exponent
    return n0 + n1 * abs(i)

def sigma(theta, c, k_c_prime, rho, g, b, k_phi_prime, r, n, theta_1, theta_m, theta_2): #Normal Stress on the Wheel
    common_factor = (c * k_c_prime + rho * g * b * k_phi_prime) * (r / b) ** n
    if theta_m < theta < theta_1:
        return common_factor * (np.cos(theta) - np.cos(theta_1)) ** n
    elif theta_2 < theta < theta_m:
        return common_factor * (np.cos(theta_1 - (theta - theta_2) / (theta_m - theta_2) * (theta_1 - theta_m)) - np.cos(theta_1)) ** n
    else:
        return 0

def j(r_s, theta_1, theta, i): #Shear Displacement
    if i >= 0:
        return r_s * (theta_1 - theta - (1 - i) * (np.sin(theta_1) - np.sin(theta)))
    else:
        return r_s * (theta_1 - theta - (np.sin(theta_1) - np.sin(theta)) / (1 + i))

def tau_max(c, theta, phi, sigma_value): #Maximum Shear Stress on the Wheel
    return c + sigma_value * np.tan(phi)

def tau(theta, tau_max_value, j_value, k): #Shear Stress on the Wheel
    return tau_max_value * (1 - np.exp(-j_value / k))

def compute_w(params): #Vertical Load on the Wheel
    p = params.get_all_params()
    p['r_s'] = calculate_r_s(p['r'], p['h'])

    i_value = calculate_i(p['v'], p['r_s'], p['omega'])
    n_value = calculate_n(p['n0'], p['n1'], i_value)
    theta_m_value = calculate_theta_m(p['a1'], p['a2'], i_value, p['theta_1'])

    try:
        # ---  integrate σ(θ) · cosθ ---
        integral_sigma_cos = quad(lambda theta: sigma(theta, p['c'], p['k_c_prime'], p['rho'], p['g'],p['b'], p['k_phi_prime'], p['r'],n_value, p['theta_1'], theta_m_value, p['theta_2']) * np.cos(th),p['theta_2'], p['theta_1'])[0]

        # --- integrate τ(θ) · sinθ  ---
        integral_tau_sin = quad(lambda theta: tau(theta,tau_max(p['c'], th, p['phi'],sigma(theta, p['c'], p['k_c_prime'], p['rho'], p['g'],p['b'], p['k_phi_prime'], p['r'],n_value, p['theta_1'], theta_m_value, p['theta_2'])),j(p['r_s'], p['theta_1'], th, i_value),p['k']) * np.sin(theta),p['theta_2'], p['theta_1'])[0]

        # --- Assemble F_W ---
        W = p['b'] * (p['r'] * integral_sigma_cos + p['r_s'] * integral_tau_sin)
        return W
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

def compute_d(params): #Drawbar Pull Force
    p = params.get_all_params()
    p['r_s'] = calculate_r_s(p['r'], p['h'])

    i_value = calculate_i(p['v'], p['r_s'], p['omega'])
    n_value = calculate_n(p['n0'], p['n1'], i_value)
    theta_m_value = calculate_theta_m(p['a1'], p['a2'], i_value, p['theta_1'])

    try:
        # integrate τ(θ)·cosθ
        int_tau_cos = quad(lambda theta: tau(theta,tau_max(p['c'], theta, p['phi'],sigma(th, p['c'], p['k_c_prime'], p['rho'], p['g'],p['b'], p['k_phi_prime'], p['r'],n_value, p['theta1'], theta_m_value, p['theta2'])),j(p['r_s'], p['theta1'], theta, i_value), p['k']) * np.cos(theta),p['theta2'], p['theta1'])[0]

        # integrate σ(θ)·sinθ
        int_sigma_sin = quad(lambda theta: sigma(theta, p['c'], p['k_c_prime'], p['rho'], p['g'],p['b'], p['k_phi_prime'], p['r'],n_value, p['theta1'], theta_m_value, p['theta2']) * np.sin(theta),p['theta2'], p['theta1'])[0]

        F_DP = p['b'] * ( p['r_s'] * int_tau_cos - p['r'] * int_sigma_sin )
        return D
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

def compute_t(params): #Resistive Moment on the Wheel
    p = params.get_all_params()
    p['r_s'] = calculate_r_s(p['r'], p['h'])

    i_value = calculate_i(p['v'], p['r_s'], p['omega'])
    n_value = calculate_n(p['n0'], p['n1'], i_value)
    theta_m_value = calculate_theta_m(p['a1'], p['a2'], i_value, p['theta_1'])

    try:
        integral2 = quad(lambda theta: tau(theta, tau_max(p['c'], theta, p['phi'], sigma(theta, p['c'], p['k_c_prime'], p['rho'], p['g'], p['b'], p['k_phi_prime'], p['r'], n_value, p['theta_1'], theta_m_value, p['theta_2'])),
                                           j(p['r_s'], p['theta_1'], theta, i_value), p['k']), p['theta_2'], p['theta_1'])[0]

        T = (p['r_s'] ** 2)* p['b'] * (integral2)
        return T
    except Exception as e:
        print(f"An error occurred: {e}")
        return None
