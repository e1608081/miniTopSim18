# -*- coding: utf-8 -*-
"""
Surface object and functions for miniTopSim
"""


import numpy as np
import matplotlib.pyplot as plt

from init_surface import init_surface


class Surface:
    
    def __init__(self):
        '''Init the object. Copy start value for printing.'''
        self.x = np.linspace(-50, 50, 101)
        self.y = np.array(init_surface(self.x))
        self.startx = np.copy(self.x)
        self.starty = np.copy(self.y)
        
    def normal(self):
        '''Calc normal vector and return it.'''        
        x = self.enlarge(self.x)
        y = self.enlarge(self.y)
        dx = self.calc_vector(self.x, x)
        dy = self.calc_vector(self.y, y)
        length = np.linalg.norm([dx,dy],axis=0)
        return dy/length, -dx/length
    
    def enlarge(self, value):
        '''Enlarge array to calc nomal vector.'''
        new_value = np.zeros(value.size + 2)
        new_value[1:-1] = value
        new_value[0] = new_value[1]
        new_value[-1] = new_value[-2]
        return new_value
    
    def calc_vector(self, orig, value):
        '''Subtract end coordinate with start coordinate.'''
        delta = value[2:] - value[:-2]
        return delta
        
    def plot(self, tend, dtime):
        '''Plot figure with start and end values.'''
        fig = plt.figure()
        plt.title('miniTopSim')
        plt.grid()
        plt.plot(self.startx, self.starty, 'gx-', label='Start')
        plt.plot(self.x, self.y, 'rx-', label='Stop')
        plt.legend()
        plt.show()
        fig.savefig('basic_{}_{}.png'.format(int(tend), int(dtime)), dpi=500)
        
    def write(self, file, time):  
        '''Write output file.'''
        with open(file, 'a') as fp:
            fp.write('surface: {}, {}, x-positions y-positions\n'.format(time, len(self.x)))
            for x, y in zip(self.x, self.y):
                fp.write("{} {}\n".format(x, y))