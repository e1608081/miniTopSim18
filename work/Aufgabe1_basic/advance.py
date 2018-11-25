# -*- coding: utf-8 -*-
"""
Functions to calculate time and coordinates
"""


def advance(surface, dtime):    
    '''Calculate coordinates after each step.'''           
    normal_x, normal_y = surface.normal()
    surface.x +=  normal_x*dtime
    surface.y +=  normal_y*dtime
     
    
def timestep(dtime, time, endTime):
    '''Get next possible timestep.'''
    if time + dtime < endTime:
        return dtime     
    else:
        return endTime - time