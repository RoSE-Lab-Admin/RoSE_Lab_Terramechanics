# Update RoverSimulator to support individual wheel running states from their own parameter classes

import numpy as np
from .Rover_Class import Rover
from .Wheel_Class import Wheel
from .Parameter_Manager_Split import RunningStateParameters

class RoverSimulator:
    def __init__(self, mass, length, width, height, soil_params, wheel_params):
        self.mass = mass  # Total rover mass (kg)
        self.length = length  # Distance front to back (m)
        self.width = width    # Distance left to right (m)
        self.height = height  # Height of center of mass (m)
        self.gravity = 9.8    # m/s^2

        self.rover = Rover(soil_params, wheel_params, wheel_count=4)

        # Initial kinematic state
        self.position = np.array([0.0, 0.0])  # x, y
        self.velocity = np.array([0.0, 0.0])  # vx, vy

        self.slope_angle = 0.0  # radians

    def set_wheel_inputs(self, omega_list, theta_1_list, theta_2_list):
        """Set individual wheel velocities and angles directly via RunningStateParameters."""
        self.rover.update_wheel_states(omega_list, theta_1_list, theta_2_list)

    def set_slope_angle(self, angle_rad):
        """Set terrain slope angle (in radians) relative to the global x-axis."""
        self.slope_angle = angle_rad

    def simulate_step(self, dt):
        """Simulates one time-step of rover motion."""
        # Get forces from all wheels
        forces = self.rover.get_total_forces()
        total_drawbar_pull = forces['Total_F_DP']

        # Gravity force in body frame (accounts for slope)
        gravity_vector = np.array([0, -self.mass * self.gravity])
        rotation_matrix = np.array([
            [np.cos(self.slope_angle), -np.sin(self.slope_angle)],
            [np.sin(self.slope_angle),  np.cos(self.slope_angle)]
        ])
        gravity_global = rotation_matrix @ gravity_vector

        # Net force = drawbar pull (x-direction) + gravity
        net_force = np.array([total_drawbar_pull, 0.0]) + gravity_global

        # Acceleration and velocity update
        acceleration = net_force / self.mass
        self.velocity += acceleration * dt
        self.position += self.velocity * dt

        return {
            'position': self.position.copy(),
            'velocity': self.velocity.copy(),
            'acceleration': acceleration.copy(),
            'drawbar_pull': total_drawbar_pull,
            'gravity_force': gravity_global.copy(),
            'net_force': net_force.copy()
        }
