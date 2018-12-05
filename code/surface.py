# -*- coding: utf-8 -*-
"""
Surface object and functions for miniTopSim

Class functions:
    normal(): calculates the normal vector and returns it
    plot(tend, dtime): plots the start and given values
    write(file, time): writes output file
"""


import numpy as np
import matplotlib.pyplot as plt

from init_surface import init_surface
import parameters as par

class Surface:
    
    def __init__(self):
        """Init the object. Copy start value for printing."""
        num_points = (par.XMAX - par.XMIN) // par.DELTA_X
        self.x = np.linspace(par.XMIN, par.XMAX, num_points)
        self.y = np.array(init_surface(self.x))
        self.startx = np.copy(self.x)
        self.starty = np.copy(self.y)
        
    def normal(self):
        """Calc normal vector and return it."""        
        x = np.concatenate(([self.x[0]], self.x, [self.x[-1]]))
        y = np.concatenate(([self.y[0]], self.y, [self.y[-1]]))
        dx = self._calc_vector(x)
        dy = self._calc_vector(y)
        length = np.linalg.norm([dx,dy],axis=0)
        return dy/length, -dx/length
    
    def _calc_vector(self, value):
        """Subtract end coordinate with start coordinate.
        
        :param value: coordinate array
        """
        delta = value[2:] - value[:-2]
        return delta
               
    def write(self, file, time, mode):  
        """Write output file.
        
        :param file: filename
        :param time: current time
        """
        with open(file, mode) as fp:
            fp.write('surface: {} {} x-positions y-positions\n'.format(time, len(self.x)))
            for x, y in zip(self.x, self.y):
                fp.write("{} {}\n".format(x, y))