class ParameterManager:
    def __init__(self):
        self.params = {
            'v': {'value': 0.057, 'bounds': (0, 0.14), 'uncertainty': 0.024, 'name': 'Wheel Linear Velocity', 'unit': 'm/s'},
            'r': {'value': 0.104, 'bounds': (0, 0.25), 'uncertainty': 0.027, 'name': 'Wheel Radius', 'unit': 'm'},
            'omega': {'value': 0.699, 'bounds': (0.1, 2), 'uncertainty': 0.401, 'name': 'Wheel Angular Velocity', 'unit': 'rad/s'},
            'a1': {'value': 0.42, 'bounds': (0.0, 1), 'uncertainty': 0.07, 'name': 'Model Parameter a1', 'unit': 'unitless'},
            'a2': {'value': 0.38, 'bounds': (0.0, 1), 'uncertainty': 0.2, 'name': 'Model Parameter a2', 'unit': 'unitless'},
            'n0': {'value': 1.23, 'bounds': (1.0, 2.0), 'uncertainty': 0.24, 'name': 'Nominal Sinkage Exponent', 'unit': 'unitless'},
            'n1': {'value': 0.0045, 'bounds': (-1, 1), 'uncertainty': 0.49, 'name': 'Slip-Sinkage Exponent', 'unit': 'unitless'},
            'b': {'value': 0.08, 'bounds': (0, 1), 'uncertainty': 0.043, 'name': 'Wheel Width', 'unit': 'm'},
            'theta_1': {'value': 0.53860861, 'bounds': (0, 1.5708), 'uncertainty': 0.42289328, 'name': 'Wheel Entry Angle', 'unit': 'rad'},
            'theta_2': {'value': -0.08464847, 'bounds': (-1.5708, 0), 'uncertainty': 0.1197296, 'name': 'Wheel Exit Angle', 'unit': 'rad'},
            'c': {'value': 0.68, 'bounds': (0, 2), 'uncertainty': 0.29, 'name': 'Cohesion', 'unit': 'kPa'},
            'k_c_prime': {'value': 100.41, 'bounds': (0, 1000), 'uncertainty': 254.53, 'name': 'Bekker Coefficient of Cohesion\'', 'unit': 'kPa/m^n'},
            'rho': {'value': 1257.65, 'bounds': (0, 3000), 'uncertainty': 849.78, 'name': 'Soil Mass Density', 'unit': 'kg/m^3'},
            'g': {'value': 9.8, 'bounds': (9.8, 9.8), 'uncertainty': 0, 'name': 'Gravitational Acceleration', 'unit': 'm/s^2'},  # Constant
            'k_phi_prime': {'value': 945.85, 'bounds': (0, 3000), 'uncertainty': 903.85, 'name': 'Bekker Coefficient of Density\'', 'unit': 'kPa/m^n-1'},
            'phi': {'value': 0.59969513, 'bounds': (0, 1.5708), 'uncertainty': 0.08272861, 'name': 'Internal Friction Angle', 'unit': 'rad'},
            'k': {'value': 0.2, 'bounds': (0.0, 1), 'uncertainty': 0.4, 'name': 'Shear Deformation Modulus', 'unit': 'm'},
            'h': {'value': 0.008, 'bounds': (0.001, 0.004), 'uncertainty': 0.002, 'name': 'Wheel Grouser height', 'unit': 'm'}
        }

        # Update bounds based on uncertainty, clamping to physical bounds
        for param, details in self.params.items():
            uncertainty = details['uncertainty']
            value = details['value']
            phys_lower_bound, phys_upper_bound = details['bounds']
            lower_bound = max(phys_lower_bound, value - 10 * uncertainty)
            upper_bound = min(phys_upper_bound, value + 10 * uncertainty)
            details['bounds'] = (lower_bound, upper_bound)

    def get(self, param_name):
        return self.params[param_name]['value']

    def set(self, param_name, value):
        self.params[param_name]['value'] = value

    def get_bounds(self, param_name):
        return self.params[param_name]['bounds']

    def get_uncertainty(self, param_name):
        return self.params[param_name]['uncertainty']

    def get_name(self, param_name):
        return self.params[param_name]['name']

    def get_unit(self, param_name):
        return self.params[param_name]['unit']

    def get_all_params(self):
        return {key: self.params[key]['value'] for key in self.params.keys()}

# Initialize parameter manager
params = ParameterManager()
