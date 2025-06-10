class Wheel:
    def __init__(self, soil_params, wheel_params, running_state):
        self.soil_params = soil_params
        self.wheel_params = wheel_params
        self.running_state = running_state

    def get_combined_params(self):
        params = {}
        params.update(self.soil_params.get_all_params())
        params.update(self.wheel_params.get_all_params())
        params.update(self.running_state.get_all_params())
        return params

    def compute_forces(self):
        class Wrapper:
            def __init__(self, full_params):
                self.full_params = full_params

            def get_all_params(self):
                return self.full_params

        wrapped_params = Wrapper(self.get_combined_params())
        F_W = compute_w(wrapped_params)
        F_DP = compute_d(wrapped_params)
        T = compute_t(wrapped_params)
        return {
            'F_W': F_W,
            'F_DP': F_DP,
            'T': T
        }

    def update_running_state(self, omega, theta_1, theta_2):
        self.running_state.params['omega']['value'] = omega
        self.running_state.params['theta_1']['value'] = theta_1
        self.running_state.params['theta_2']['value'] = theta_2
