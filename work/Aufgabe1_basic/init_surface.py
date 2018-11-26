# -*- coding: utf-8 -*-
"""
Init function for surface
"""


import numpy as np


def init_surface(x):
    '''calculate start coordinates after formula. Crop values at +/- 25.'''
    mask = np.where(np.abs(x) < 25)
    y = np.zeros(x.size)
    y[mask] = -50 * (1 + np.cos((2 * np.pi * x[mask]) / 50))
    return y