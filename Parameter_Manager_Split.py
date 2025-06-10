# Refactored ParameterManager with no uncertainties and grouped by Ps, Pw, Pr

class SoilParameters:
    def __init__(self):
        self.params = {
            'a1': {'value': 0.42, 'bounds': (0.0, 1), 'name': 'Model Parameter a1', 'unit': 'unitless'},
            'a2': {'value': 0.38, 'bounds': (0.0, 1), 'name': 'Model Parameter a2', 'unit': 'unitless'},
            'n0': {'value': 1.23, 'bounds': (1.0, 2.0), 'name': 'Nominal Sinkage Exponent', 'unit': 'unitless'},
            'n1': {'value': 0.0045, 'bounds': (-1, 1), 'name': 'Slip-Sinkage Exponent', 'unit': 'unitless'},
            'c': {'value': 0.68, 'bounds': (0, 2), 'name': 'Cohesion', 'unit': 'kPa'},
            'k_c_prime': {'value': 100.41, 'bounds': (0, 1000), 'name': "Bekker Coefficient of Cohesion'", 'unit': 'kPa/m^n'},
            'rho': {'value': 1257.65, 'bounds': (0, 3000), 'name': 'Soil Mass Density', 'unit': 'kg/m^3'},
            'k_phi_prime': {'value': 945.85, 'bounds': (0, 3000), 'name': "Bekker Coefficient of Density'", 'unit': 'kPa/m^n-1'},
            'phi': {'value': 0.59969513, 'bounds': (0, 1.5708), 'name': 'Internal Friction Angle', 'unit': 'rad'},
            'k': {'value': 0.2, 'bounds': (0.0, 1), 'name': 'Shear Deformation Modulus', 'unit': 'm'},
        }

    def get_all_params(self):
        return {k: v['value'] for k, v in self.params.items()}


class WheelParameters:
    def __init__(self):
        self.params = {
            'v': {'value': 0.057, 'bounds': (0, 0.14), 'name': 'Wheel Linear Velocity', 'unit': 'm/s'},
            'r': {'value': 0.104, 'bounds': (0, 0.25), 'name': 'Wheel Radius', 'unit': 'm'},
            'b': {'value': 0.08, 'bounds': (0, 1), 'name': 'Wheel Width', 'unit': 'm'},
            'g': {'value': 9.8, 'bounds': (9.8, 9.8), 'name': 'Gravitational Acceleration', 'unit': 'm/s^2'},
            'h': {'value': 0.008, 'bounds': (0.001, 0.004), 'name': 'Wheel Grouser Height', 'unit': 'm'},
        }

    def get_all_params(self):
        return {k: v['value'] for k, v in self.params.items()}


class RunningStateParameters:
    def __init__(self):
        self.params = {
            'omega': {'value': 0.699, 'bounds': (0.1, 2), 'name': 'Wheel Angular Velocity', 'unit': 'rad/s'},
            'theta_1': {'value': 0.53860861, 'bounds': (0, 1.5708), 'name': 'Wheel Entry Angle', 'unit': 'rad'},
            'theta_2': {'value': -0.08464847, 'bounds': (-1.5708, 0), 'name': 'Wheel Exit Angle', 'unit': 'rad'},
        }

    def get_all_params(self):
        return {k: v['value'] for k, v in self.params.items()}


# Combined Parameter Manager for full access
class ParameterManager:
    def __init__(self):
        self.soil = SoilParameters()
        self.wheel = WheelParameters()
        self.running_state = RunningStateParameters()

    def get_all_params(self):
        all_params = {}
        all_params.update(self.soil.get_all_params())
        all_params.update(self.wheel.get_all_params())
        all_params.update(self.running_state.get_all_params())
        return all_params

# Initialize parameter manager
params = ParameterManager()
