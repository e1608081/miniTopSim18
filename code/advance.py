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

import parameters as par
from sputtering import sputter_yield
import surface 


def advance(surface, dtime):
    """Calculate coordinates after each step.

    :param surface: surface object
    :param dtime: timestep size
    """

    normal_v = get_velocities(surface, dtime)
    normal_x, normal_y = surface.normal()

    if par.TIME_INTEGRATION == 'normal':
        surface.x += dtime * normal_x * normal_v
        surface.y += dtime * normal_y * normal_v

    if par.TIME_INTEGRATION == 'vertical':
        surface.y += dtime * normal_v / normal_y
    
    surface.deloop()
    
    surface.eliminate_overhangs()

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
    """ Calculates the surface velocity depending on [SETUP] ETCHING.

    In case of sputtering the angle between the beam direction and surface
    normal is calculated and then the sputter_yield() function is called
    from the sputterung-module.

    :return: surface velocity [unit: nm/s]
    """

    if par.ETCHING is True:
        # If etching is used
        return par.ETCH_RATE * np.ones_like(surface.x)
    else:
        # beam current density J = F_beam * e
        J = par.BEAM_CURRENT_DENSITY
        N = par.DENSITY

        F_beam = J / e

        normal_x, normal_y = surface.normal()

        theta = np.arccos(-np.minimum(normal_y, 0.))

        Y_s = sputter_yield(theta)

        # F_sput = F_beam * Y_s(theta) * cos(theta)
        F_sput = F_beam * Y_s * np.cos(theta)

        if False:
            # remove overhanging structures
            lastx = surface.x[0]
            for i in range(1, surface.x.size-1):
                if surface.x[i] < 0. and surface.x[i] < lastx:
                    if surface.y[i] <= surface.y[i-1]:
                        # shadowed point
                        F_sput[i] = 0
                    else:
                        # point that shadow the last point
                        if surface.y[i] > np.max(surface.y):
                            # update lastx to the current point
                            lastx = surface.x[i]
                        F_sput[i-1] = 0
                else:
                    lastx = surface.x[i]

            lastx = surface.x[-1]
            for i in range(-2, -surface.x.size+1, -1):
                if surface.x[i] > 0. and surface.x[i] > lastx:
                    if surface.y[i] <= surface.y[i+1]:
                        # shadowed point
                        F_sput[i] = 0
                    else:
                        # point that shadow the last point
                        if surface.y[i] > np.max(surface.y):
                            # update lastx to the current point
                            lastx = surface.x[i]
                        F_sput[i+1] = 0
                else:
                    lastx = surface.x[i]

        v_normal = F_sput / N

        # convert units from cm/s to nm/s
        v_normal = v_normal * 1e7

        return v_normal

