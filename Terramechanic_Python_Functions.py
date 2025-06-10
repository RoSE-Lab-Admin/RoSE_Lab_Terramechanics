# Updated and corrected version of the terramechanics model functions
import numpy as np
from scipy.integrate import quad

# Effective Shearing Radius
def calculate_r_s(r, h):  
    return r + 0.5 * h  # Assuming λ_s = 0.5

# Wheel Slip Ratio
def calculate_i(v, r_s, omega):  
    if v < r_s * omega:
        return 1 - (v / (r_s * omega))
    elif v > r_s * omega:
        return (r_s * omega / v) - 1
    else:
        return 0

# Angle of Maximum Normal Stress
def calculate_theta_m(a1, a2, i, theta_1):  
    return (a1 + a2 * i) * theta_1

# Sinkage Exponent based on slip ratio
def calculate_n(n0, n1, i):  
    return n0 + n1 * abs(i)

# Normal Stress Distribution
def sigma(theta, c, k_c_prime, rho, g, b, k_phi_prime, r, n, theta_1, theta_m, theta_2):  
    common_factor = (c * k_c_prime + rho * g * b * k_phi_prime) * (r / b) ** n
    if theta_m < theta < theta_1:
        return common_factor * (np.cos(theta) - np.cos(theta_1)) ** n
    elif theta_2 < theta < theta_m:
        inner_angle = theta_1 - ((theta - theta_2) / (theta_m - theta_2)) * (theta_1 - theta_m)
        return common_factor * (np.cos(inner_angle) - np.cos(theta_1)) ** n
    else:
        return 0

# Shear Displacement along the wheel-soil interface
def j(r_s, theta_1, theta, i):  
    if i >= 0:
        return r_s * (theta_1 - theta - (1 - i) * (np.sin(theta_1) - np.sin(theta)))
    else:
        return r_s * (theta_1 - theta - (np.sin(theta_1) - np.sin(theta)) / (1 + i))

# Maximum Shear Stress
def tau_max(c, theta, phi, sigma_value):  
    return c + sigma_value * np.tan(phi)

# Shear Stress as a function of shear displacement
def tau(theta, tau_max_value, j_value, k):  
    return tau_max_value * (1 - np.exp(-j_value / k))

# Vertical Load on the Wheel
def compute_vertical_load(params):  
    p = params.get_all_params()
    p['r_s'] = calculate_r_s(p['r'], p['h'])

    i_value = calculate_i(p['v'], p['r_s'], p['omega'])
    n_value = calculate_n(p['n0'], p['n1'], i_value)
    theta_m_value = calculate_theta_m(p['a1'], p['a2'], i_value, p['theta_1'])

    try:
        integral1 = quad(lambda theta: sigma(theta, p['c'], p['k_c_prime'], p['rho'], p['g'], p['b'], p['k_phi_prime'], p['r'], n_value, p['theta_1'], theta_m_value, p['theta_2']) * np.cos(theta),
                         p['theta_2'], p['theta_1'])[0]
        integral2 = quad(lambda theta: tau(theta, tau_max(p['c'], theta, p['phi'], sigma(theta, p['c'], p['k_c_prime'], p['rho'], p['g'], p['b'], p['k_phi_prime'], p['r'], n_value, p['theta_1'], theta_m_value, p['theta_2'])),
                                           j(p['r_s'], p['theta_1'], theta, i_value), p['k']) * np.sin(theta),
                         p['theta_2'], p['theta_1'])[0]

        W = p['b'] * (p['r'] * integral1 + p['r_s'] * integral2)
        return W
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

# Drawbar Pull Force
def compute_drawbar_pull(params):  
    p = params.get_all_params()
    p['r_s'] = calculate_r_s(p['r'], p['h'])

    i_value = calculate_i(p['v'], p['r_s'], p['omega'])
    n_value = calculate_n(p['n0'], p['n1'], i_value)
    theta_m_value = calculate_theta_m(p['a1'], p['a2'], i_value, p['theta_1'])

    try:
        integral1 = quad(lambda theta: sigma(theta, p['c'], p['k_c_prime'], p['rho'], p['g'], p['b'], p['k_phi_prime'], p['r'], n_value, p['theta_1'], theta_m_value, p['theta_2']) * np.sin(theta),
                         p['theta_2'], p['theta_1'])[0]
        integral2 = quad(lambda theta: tau(theta, tau_max(p['c'], theta, p['phi'], sigma(theta, p['c'], p['k_c_prime'], p['rho'], p['g'], p['b'], p['k_phi_prime'], p['r'], n_value, p['theta_1'], theta_m_value, p['theta_2'])),
                                           j(p['r_s'], p['theta_1'], theta, i_value), p['k']) * np.cos(theta),
                         p['theta_2'], p['theta_1'])[0]

        D = p['b'] * (p['r_s'] * integral2 - p['r'] * integral1)
        return D
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

# Resistive Moment on the Wheel
def compute_resustive_moment(params):  
    p = params.get_all_params()
    p['r_s'] = calculate_r_s(p['r'], p['h'])

    i_value = calculate_i(p['v'], p['r_s'], p['omega'])
    n_value = calculate_n(p['n0'], p['n1'], i_value)
    theta_m_value = calculate_theta_m(p['a1'], p['a2'], i_value, p['theta_1'])

    try:
        integral = quad(lambda theta: tau(theta, tau_max(p['c'], theta, p['phi'], sigma(theta, p['c'], p['k_c_prime'], p['rho'], p['g'], p['b'], p['k_phi_prime'], p['r'], n_value, p['theta_1'], theta_m_value, p['theta_2'])),
                                         j(p['r_s'], p['theta_1'], theta, i_value), p['k']), 
                        p['theta_2'], p['theta_1'])[0]

        T = p['b'] * (p['r_s'] ** 2) * integral
        return T
    except Exception as e:
        print(f"An error occurred: {e}")
        return None
