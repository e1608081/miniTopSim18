# -*- coding: utf-8 -*-

"""
Functions to calculate time and coordinates

Module functions:
    advance(surface, dtime)
    timestep(dtime, time, endTime)
    get_velocities(surface, dtime)
"""

import numpy as np
from scipy.constants import e
from scipy.spatial.distance import cdist

import parameters as par
from sputtering import sputter_yield, sputter_yield_derivative
import surface
from beam import Beam


def advance(surface, dtime):
    """Calculate coordinates after each step.

    Use either normal or vertical calculation, according to the config.
    'normal' moves the surface points in direction of the normal vector,
    'vertical' uses only the y direction and leaves the x coordinate
    fixed.

    :param surface: surface object
    :param dtime: timestep size
    """

    normal_v, d_normal_v = get_velocities(surface, dtime)
    normal_x, normal_y = surface.normal()

    if par.TIME_INTEGRATION == 'normal':
        surface.x += dtime * normal_x * normal_v
        surface.y += dtime * normal_y * normal_v

    elif par.TIME_INTEGRATION == 'vertical':
        surface.y += dtime * normal_v / normal_y

    elif par.TIME_INTEGRATION == 'characteristics':
        theta = np.arccos(-np.minimum(normal_y, 0.))
        #theta is positiv if the slope is positiv
        theta = np.copysign(theta, normal_x)

        v_x = normal_v * np.sin(theta) + d_normal_v * np.cos(theta)
        v_y = -normal_v * np.cos(theta) + d_normal_v * np.sin(theta)

        surface.x += v_x * dtime
        surface.y += v_y * dtime

    surface.deloop()
    surface.eliminate_overhangs()

    if par.ADAPTIVE_GRID is True:
        surface.adapt()

def timestep(dtime, time, endTime):
    """Get next possible timestep.

    :param dtime: timestep size
    :param time: current time
    :param endTime: end time
    """
    if time + dtime < endTime:
        return dtime
    else:
        return endTime - time


def get_velocities(surface, dtime):
    """ Calculates the surface velocity and derivatve depending on [SETUP] ETCHING.

    In case of sputtering the angle between the beam direction and surface
    normal is calculated and then the sputter_yield() function is called
    from the sputterung-module.

    :return: tuple of surface velocity [unit: nm/s] and the derivative
    """

    if par.ETCHING is True:
        # If etching is used
        return par.ETCH_RATE * np.ones_like(surface.x), 0
    else:
        N = par.DENSITY

        beam = Beam()
        F_beam = beam(surface.x)

        normal_x, normal_y = surface.normal()

        theta = np.arccos(-np.minimum(normal_y, 0.))
        #theta is positiv if the slope is positiv
        theta = np.copysign(theta, normal_x)

        Y_s = sputter_yield(theta)

        # F_sput = F_beam * Y_s(theta) * cos(theta)
        F_sput = F_beam * Y_s * np.cos(theta)


        if par.REDEP is True:
            viewfactor = surface.calc_viewfactor()
            F_redep = viewfactor.dot(F_sput)
        else:
            F_redep = 0

        v_normal = (F_sput - F_redep) / N

        # convert units from cm/s to nm/s
        v_normal = v_normal * 1e7


        #claculate derivative
        dY_s = sputter_yield_derivative(theta)
        dF_sput = F_beam * (dY_s * np.cos(theta) - Y_s * np.sin(theta))

        if par.REDEP is True:
            viewfactor = surface.calc_viewfactor_derivative()
            dF_redep = viewfactor.dot(F_sput)
        else:
            dF_redep = 0

        d_v_normal = (dF_sput - dF_redep) / N
        # convert units from cm/s^2 to nm/s^2
        d_v_normal = d_v_normal * 1e7


        return v_normal, d_v_normal
