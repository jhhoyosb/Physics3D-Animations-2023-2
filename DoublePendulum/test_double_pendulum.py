"""
Test of the double pendulum using normal nodes
"""
from double_pendulum import rk4
import numpy as np


def test_double_pendulum():
    """
    test of rk4
    :return:
    """
    # define constants
    gravity = 9.8
    length = 0.5
    wsum = np.sqrt(2 + np.sqrt(2)) * np.sqrt(gravity / length)
    wminus = np.sqrt(2 - np.sqrt(2)) * np.sqrt(gravity / length)
    array_time = np.arange(0, 5 + 0.005, 0.005)
    constant_a = np.deg2rad(5)
    constant_b = 0

    # Calculo de los thetas
    theta1 = constant_a * np.cos(wsum * array_time) + constant_b * np.cos(wminus * array_time)
    theta2 = (constant_a * -np.sqrt(2) * np.cos(wsum * array_time) +
              constant_b * np.sqrt(2) * np.cos(wminus * array_time))

    # Numeric method results
    initial_condition_1 = np.array([np.deg2rad(5), 0])
    initial_condition_2 = np.array([-np.deg2rad(5 * np.sqrt(2)), 0])

    array_time, numeric_theta1, numeric_theta2 = rk4(0,
                                                     5,
                                                     1000,
                                                     initial_condition_1,
                                                     initial_condition_2)

    assert abs(np.average(numeric_theta2[:, 0] - theta2)) > 1*10**-2
    assert abs(np.average(numeric_theta1[:, 0] - theta1)) < 1 * 10 ** -2
