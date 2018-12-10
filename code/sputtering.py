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


def init_sputtering():
    """ Initialize sputter yield.

    If SPUTTER_YIELD_FILE is '', then the sputter yield is calculated by the
    yamamura function otherwise by the table that is stored in ../tables/.
    """

    global get_sputter_yield

    if par.SPUTTER_YIELD_FILE == '':
        # calculate sputter yield by yamamura formula
        y0 = par.SPUTTER_YIELD_0
        f = par.SPUTTER_YIELD_F
        b = par.SPUTTER_YIELD_B
        get_sputter_yield = Sputter_yield_Yamamura(y0, f, b)
    else:
        # calculate sputter yield by table
        get_sputter_yield = Sputter_yield_table(par.SPUTTER_YIELD_FILE)


def sputter_yield(theta):
    """ Calculate the sputter yield for a theta in rad

    Important: To use this function init_sputtering() need to be run first!
    """
    return get_sputter_yield(theta)
