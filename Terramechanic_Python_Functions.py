from functools import lru_cache
from scipy.integrate import quad
import numpy as np

# ------------------------------------------------------------------
#  Helpers that change *nothing* outside their own theta dependence
# ------------------------------------------------------------------
def _sigma_cached(theta, *, c, k_c_prime, rho, g, b, k_phi_prime,
                  r, n, theta_1, theta_m, theta_2):
    """σ(θ) but with a real‑valued, safe base."""
    common = (c * k_c_prime + rho * g * b * k_phi_prime) * (r / b) ** n

    if theta_m < theta < theta_1:
        base = np.clip(np.cos(theta) - np.cos(theta_1), 0, None)
    elif theta_2 < theta < theta_m:
        alpha = (theta - theta_2) / (theta_m - theta_2)
        base = np.clip(np.cos(theta_1 - alpha * (theta_1 - theta_m)) - np.cos(theta_1),
                       0, None)
    else:
        return 0.0
    return common * base ** n


def _tau(theta, *, p_const, n_val, theta_m):
    """τ(θ) with everything except θ factored out."""
    σ = _sigma_cached(theta, **p_const, n=n_val, theta_m=theta_m)
    τ_max = p_const['c'] + σ * np.tan(p_const['phi'])
    j_val = j(p_const['r_s'], p_const['theta_1'], theta, p_const['i'])
    return τ_max * (1.0 - np.exp(-j_val / p_const['k']))


# ------------------------------------------------------------------
#  Forces & moment
# ------------------------------------------------------------------
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

def _precompute_constants(params):
    p = params.get_all_params().copy()
    p['r_s']      = calculate_r_s(p['r'], p['h'])
    p['i']        = calculate_i(p['v'], p['r_s'], p['omega'])
    p['n']        = calculate_n(p['n0'], p['n1'], p['i'])
    p['theta_m']  = calculate_theta_m(p['a1'], p['a2'], p['i'], p['theta_1'])
    return p


def compute_vertical_force(params):
    p = _precompute_constants(params)

    try:
        # ∫ σ cosθ
        f1 = lambda θ: _sigma_cached(θ, **p) * np.cos(θ)
        int1 = quad(f1, p['theta_2'], p['theta_1'])[0]

        # ∫ τ sinθ
        f2 = lambda θ: _tau(θ, p_const=p, n_val=p['n'], theta_m=p['theta_m']) * np.sin(θ)
        int2 = quad(f2, p['theta_2'], p['theta_1'])[0]

        return p['b'] * (p['r'] * int1 + p['r_s'] * int2)

    except Exception as e:
        print(f"[compute_w] {e}")
        return None


def compute_drawbar(params):
    p = _precompute_constants(params)

    try:
        f1 = lambda θ: _tau(θ, p_const=p, n_val=p['n'], theta_m=p['theta_m']) * np.cos(θ)
        int_tau_cos = quad(f1, p['theta_2'], p['theta_1'])[0]

        f2 = lambda θ: _sigma_cached(θ, **p) * np.sin(θ)
        int_sigma_sin = quad(f2, p['theta_2'], p['theta_1'])[0]

        return p['b'] * (p['r_s'] * int_tau_cos - p['r'] * int_sigma_sin)

    except Exception as e:
        print(f"[compute_d] {e}")
        return None


def compute_resistance_moment(params):
    p = _precompute_constants(params)

    try:
        f = lambda θ: _tau(θ, p_const=p, n_val=p['n'], theta_m=p['theta_m'])
        int_tau = quad(f, p['theta_2'], p['theta_1'])[0]

        return p['b'] * p['r_s'] ** 2 * int_tau

    except Exception as e:
        print(f"[compute_t] {e}")
        return None
