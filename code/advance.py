# -*- coding: utf-8 -*-
"""
Functions to calculate time and coordinates

Module functions:
    advance(surface, dtime)
    timestep(dtime, time, endTime)
"""


def advance(surface, dtime):    
    """Calculate coordinates after each step.
    
    :param surface: surface object
    :param dtime: timestep size
    """           
    normal_x, normal_y = surface.normal()
    surface.x +=  normal_x*dtime
    surface.y +=  normal_y*dtime
     
    
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