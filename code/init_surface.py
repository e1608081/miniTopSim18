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
    mask = np.where(np.logical_and(np.array(x) > par.FUN_XMIN, np.array(x) < par.FUN_XMAX))
    y = np.zeros(x.size)
    middle = (par.FUN_XMAX + par.FUN_XMIN) / 2
    len = par.FUN_XMAX - par.FUN_XMIN
    
    if par.INITIAL_SURFACE_TYPE == "Cosine":
        y[mask] = -50 * (1 + np.cos(2 * np.pi * (x[mask] - middle) / len))
    else:
        raise ValueError("use of unimplemented Function")
    
    return y