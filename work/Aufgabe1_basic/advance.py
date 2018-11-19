# -*- coding: utf-8 -*-
"""
@author: Philipp
"""


def advance(surface, dtime):               
    normal_x, normal_y = surface.normal()
    surface.x +=  normal_x*dtime
    surface.y +=  normal_y*dtime
     
    
def timestep(dtime, time, endTime):
    if time + dtime < endTime:
        return dtime     
    else:
        return endTime - time