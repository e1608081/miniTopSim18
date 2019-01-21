# -*- coding: utf-8 -*-

"""
In this module there is an initialisation function and two classes defined.
These classes can be used as callable objects.
"""

import os
import numpy as np
from scipy.interpolate import interp1d

import parameters as par


class Sputter_yield_Yamamura:
    """ Class, that calculate the sputter yield with the help of the yamamura
    function.
    """

    def __init__(self, y0, f, b):
        """ Initialisation of the parameters of the yamamura function.

        :param y0:  Parameter Y_0 of yamamura function
        :param f:   Parameter f of yamamura function
        :param b:   Parameter b of yamamura function
        """

        self.y0 = y0
        self.f = f
        self.b = b

    def __call__(self, theta):
        """ Is executed, if an object of type Sputter_yield_Yamamura is
        called.

        :param theta: in [rad]
        :return: Results of yamamura function
        """

        cos_theta = np.cos(theta)
        y = self.y0 * cos_theta**(-self.f) * np.exp(self.b - self.b/cos_theta)
        return y


class Sputter_yield_table:
    """ Class, that calculate the sputter yield with the help of a table.
    """

    def __init__(self, table_file):
        """ Read the table file and store sputter yield in an interp1d object.
        """

        filedir = os.path.dirname(__file__)
        table_file = os.path.join(filedir, '..', 'tables', table_file)
        theta, sputter_yield = np.loadtxt(table_file,
                                          skiprows=1,
                                          usecols=(0, 1, ),
                                          unpack=True)

        # self.sputter_yield_interp1d can be now interpolate for all
        # min(theta) <= theta <= max(theta).
        # For all other values of theta interp1d return a ValueError.
        # To prevent the raise of the ValueError, the bounds_error
        # attribute has to be False.
        # Now interp1d return for each other theta the fill_value attribute.

        self.sputter_yield_interp1d = interp1d(np.radians(theta),
                                               sputter_yield,
                                               bounds_error=False,
                                               fill_value=0)

    def __call__(self, theta):
        return self.sputter_yield_interp1d(theta)


class Sputter_yield_derivative:
    """Callable class that calculates the derivative of the yamamura

    In order to initalize from table file, use the get_params_from_file function
    """
    def __init__(self, y0, f, b):
        """ Initialisation of the parameters of the yamamura function.

        :param y0:  Parameter Y_0 of yamamura function
        :param f:   Parameter f of yamamura function
        :param b:   Parameter b of yamamura function
        """

        self.y0 = y0
        self.f = f
        self.b = b

    def __call__(self, theta):
        """Method where the actual derivative is calculated

        :param theta: angel in [rad]
        :return: derivative of the yamamura function dY_s/dtheta
        """

        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        exp_term = np.exp(self.b - self.b/cos_theta)

        yprime = self.y0 * self.f * cos_theta**(-self.f-1) * sin_theta * exp_term
        yprime -= self.y0 * cos_theta**(-self.f) * exp_term * self.b * sin_theta / (cos_theta**2)

        return yprime

    def get_params_from_file(table_file):
        """Approximate yamamura function

        in order to initialize the Sputter_yield_derivative class from a
        table file, the parameters first have to be approximated by the
        scipy.optimize.curve_fit method

        :params: table file
        :return: parameter for the yamamura function
        """


        from scipy import optimize

        def yamamura_opt(theta, y0, f, b):
            """ Yamamura function used for initialisation with a table
            """
            cos_theta = np.cos(theta)
            y = y0 * cos_theta**(-f) * np.exp(b - b/cos_theta)
            return y

        filedir = os.path.dirname(__file__)
        table_file = os.path.join(filedir, '..', 'tables', par.SPUTTER_YIELD_FILE)
        theta, sputter_yield = np.loadtxt(table_file,
                                          skiprows=1,
                                          usecols=(0, 1, ),
                                          unpack=True)

        params, _ = optimize.curve_fit(yamamura_opt, np.radians(theta), sputter_yield)

        return params


def init_sputtering():
    """ Initialize sputter yield and derivative.

    If SPUTTER_YIELD_FILE is '', then the sputter yield is calculated by the
    yamamura function otherwise by the table that is stored in ../tables/.
    """

    global get_sputter_yield
    global get_sputter_yield_derivative

    if par.SPUTTER_YIELD_FILE == '':
        # calculate sputter yield by yamamura formula
        y0 = par.SPUTTER_YIELD_0
        f = par.SPUTTER_YIELD_F
        b = par.SPUTTER_YIELD_B
        get_sputter_yield = Sputter_yield_Yamamura(y0, f, b)
        get_sputter_yield_derivative = Sputter_yield_derivative(y0, f, b)
    else:
        # calculate sputter yield by table
        get_sputter_yield = Sputter_yield_table(par.SPUTTER_YIELD_FILE)

        # initialize sputteryield derivative from table file
        params = Sputter_yield_derivative.get_params_from_file(par.SPUTTER_YIELD_FILE)
        get_sputter_yield_derivative = Sputter_yield_derivative(*params)


def sputter_yield(theta):
    """ Calculate the sputter yield for a theta in rad

    Important: To use this function init_sputtering() need to be run first!
    """
    return get_sputter_yield(theta)


def sputter_yield_derivative(theta):
    """ Calculate the derivative of the sputter yield for theta in radians

    Important: To use this function init_sputtering() need to be run first!
    """
    return get_sputter_yield_derivative(theta)
