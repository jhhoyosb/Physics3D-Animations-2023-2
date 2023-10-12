"""
double pendulum calculation
"""
# First import libraries
import numpy as np
import matplotlib.pyplot as plt


# equations of theta dot
def equations(initial_conditions1, initial_conditions_2):
    """
    :param initial_conditions1:
    :param initial_conditions_2:
    :return: time, and values for theta 1 and 2 , and its derivatives
    """
    size = len(initial_conditions1)
    array_theta_2 = np.zeros(size)
    array_theta_1 = np.zeros(size)
    # constants
    gravity = 9.8
    mass_1 = 1
    mass_2 = 1
    length_1 = 0.5
    length_2 = 0.5
    array_theta_1[0] = initial_conditions1[1]
    array_theta_2[0] = initial_conditions_2[1]

    array_theta_1[1] = -((gravity * (2 * mass_1 + mass_2) * np.sin(initial_conditions1[0]) +
                          mass_2 * (gravity * np.sin(initial_conditions1[0] -
                                                     2 * initial_conditions_2[0]) + 2 * (
                            length_2 * initial_conditions_2[1] ** 2 +
                            length_1 * initial_conditions1[1] ** 2 * np.cos(
                                initial_conditions1[0] -
                                initial_conditions_2[0])) * np.sin(
                            initial_conditions1[0] - initial_conditions_2[0]))) / (
                                 2 * length_1 * (mass_1 + mass_2 - mass_2 * np.cos(initial_conditions1[0] -
                                                                                   initial_conditions_2[0]) ** 2)))

    array_theta_2[1] = (((mass_1 + mass_2) * (length_1 * initial_conditions1[1] ** 2 +
                                              gravity * np.cos(initial_conditions1[0])) + length_2 * mass_2 *
                         initial_conditions_2[1] ** 2 * np.cos(initial_conditions1[0] -
                        initial_conditions_2[0])) *
                        np.sin(initial_conditions1[0] - initial_conditions_2[0])) / (
                               length_2 * (mass_1 + mass_2 - mass_2 * np.cos(initial_conditions1[0] -
                                                                             initial_conditions_2[0]) ** 2))
    return array_theta_1, array_theta_2


# rk4 numeric method used to calculate both theta's
def rk4(initial_time, end_time, steps, initial_conditions1, initial_conditions2):
    """
    :param initial_time: 0 normally
    :param end_time: time in seconds
    :param steps: how many steps
    :param initial_conditions1:
    :param initial_conditions2:
    :return: values for the time array and theta 1 and 2 , and its derivatives
    """
    delta_time = (end_time - initial_time) / steps
    time_array = np.arange(initial_time, end_time + delta_time, delta_time)
    theta_1 = []
    theta_2 = []
    theta_1.append(initial_conditions1)
    theta_2.append(initial_conditions2)
    for i in range(0, steps):
        array_theta_1, array_theta_2 = equations(theta_1[i], theta_2[i])
        k1x = delta_time * array_theta_1
        k1y = delta_time * array_theta_2
        array_theta_1, array_theta_2 = equations(theta_1[i] + k1x / 2, theta_2[i] + k1y / 2)
        k2x = delta_time * array_theta_1
        k2y = delta_time * array_theta_2
        array_theta_1, array_theta_2 = equations(theta_1[i] + k2x / 2, theta_2[i] + k2y / 2)
        k3x = delta_time * array_theta_1
        k3y = delta_time * array_theta_2
        array_theta_1, array_theta_2 = equations(theta_1[i] + k3x, theta_2[i] + k3y)
        k4x = delta_time * array_theta_1
        k4y = delta_time * array_theta_2
        theta_1.append(theta_1[i] + (1 / 6) * (k1x + 2 * (k2x + k3x) + k4x))
        theta_2.append(theta_2[i] + (1 / 6) * (k1y + 2 * (k2y + k3y) + k4y))

    theta_1 = np.array(theta_1)
    theta_2 = np.array(theta_2)
    return time_array, theta_1, theta_2


u1 = np.array([np.deg2rad(5), 0])
u2 = np.array([-np.deg2rad(5 * np.sqrt(2)), 0])

time, theta1, theta2 = rk4(0, 5, 10000, u1, u2)
plt.plot(time, theta1[:, 0])
plt.plot(time, theta2[:, 0])
