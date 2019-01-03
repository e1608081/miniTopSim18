# -*- coding: utf-8 -*-
"""
Init function for surface

Module functions:
    init_surface(x): returns y value with cosine function
"""


import numpy as np
import parameters as par

def init_surface(x):
    """calculate start coordinates after formula. Crop values at +/- 25.
    
    :param x: x values of surface
    """
    middle = (par.FUN_XMAX + par.FUN_XMIN) / 2

    mask = np.where(np.logical_and(np.array(x) >= par.FUN_XMIN,
                                   np.array(x) <= par.FUN_XMAX))
    y = np.zeros_like(x)

    len = par.FUN_XMAX - par.FUN_XMIN
    
    if par.INITIAL_SURFACE_TYPE == "Cosine":
        y[mask] = par.FUN_PEAK_TO_PEAK / 2 * \
                  (1 - np.cos(2 * np.pi * (x[mask] - par.FUN_XMIN) / len))

    elif par.INITIAL_SURFACE_TYPE == 'DoubleCosine':
        y[mask] = par.FUN_PEAK_TO_PEAK / 2 * \
                  (1 - np.cos(4 * np.pi * (x[mask] - par.FUN_XMIN) / len))

    elif par.INITIAL_SURFACE_TYPE == 'Step':
        mask_left = np.where(np.array(x) < par.FUN_XMIN)
        k = -par.FUN_PEAK_TO_PEAK/len

        y[mask] = par.FUN_PEAK_TO_PEAK + k * (x[mask] - par.FUN_XMIN)
        y[mask_left] = par.FUN_PEAK_TO_PEAK

    elif par.INITIAL_SURFACE_TYPE == 'V-SHAPE':
        mask1 = np.where(np.logical_and(np.array(x) >= par.FUN_XMIN,
                                        np.array(x) <= middle))
        mask2 = np.where(np.logical_and(np.array(x) > middle,
                                        np.array(x) <= par.FUN_XMAX))
        k1 = par.FUN_PEAK_TO_PEAK / (middle - par.FUN_XMIN)
        k2 = -par.FUN_PEAK_TO_PEAK / (par.FUN_XMAX - middle)

        y[mask1] = k1 * (x[mask1] - par.FUN_XMIN)
        y[mask2] = par.FUN_PEAK_TO_PEAK + k2 * (x[mask2] - middle)

    else:
        raise ValueError("use of unimplemented Function")
    
    return y