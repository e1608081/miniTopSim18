# -*- coding: utf-8 -*-
"""
@author: Philipp
"""

import numpy as np


def init_surface(x):
    mask = np.where(np.abs(x) < 25)
    new_x = np.zeros(x.size)
    new_x[mask] = -50 * (1 + np.cos((2 * np.pi * x[mask]) / 50))
    return new_x

if __name__ == '__main__' :
    print (str(init_surface(np.linspace(-50, 50, 201))))