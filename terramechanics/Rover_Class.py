from .Parameter_Manager_Split import RunningStateParameters
from .Wheel_Class import Wheel

class Rover:
    def __init__(self, soil_params, wheel_params, wheel_count=4):
        self.wheels = []
        for _ in range(wheel_count):
            running_state = RunningStateParameters()
            wheel = Wheel(soil_params, wheel_params, running_state)
            self.wheels.append(wheel)

    def get_total_forces(self):
        total_F_W = 0
        total_F_DP = 0
        total_T = 0
        per_wheel = []
        for wheel in self.wheels:
            forces = wheel.compute_forces()
            total_F_W += forces['F_W']
            total_F_DP += forces['F_DP']
            total_T += forces['T']
            per_wheel.append(forces)
        return {
            'Total_F_W': total_F_W,
            'Total_F_DP': total_F_DP,
            'Total_T': total_T,
            'Per_Wheel': per_wheel
        }

    def update_wheel_states(self, omegas, theta_1s, theta_2s):
        for i, wheel in enumerate(self.wheels):
            wheel.update_running_state(omegas[i], theta_1s[i], theta_2s[i])
