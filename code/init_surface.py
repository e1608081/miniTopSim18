# -*- coding: utf-8 -*-
"""
Init function for surface

Module functions:
    init_surface(x): returns y value with cosine function
"""


import numpy as np


def init_surface(x):
    """calculate start coordinates after formula. Crop values at +/- 25.
    
    :param x: x values of surface
    """
    mask = np.where(np.abs(x) < 25)
    y = np.zeros(x.size)
    y[mask] = -50 * (1 + np.cos((2 * np.pi * x[mask]) / 50))
    return y