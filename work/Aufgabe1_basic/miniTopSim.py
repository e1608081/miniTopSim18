# -*- coding: utf-8 -*-
"""
Main file for miniTopSim
"""


import sys

from surface import Surface
import advance  as adv


if len(sys.argv) < 3:
    print("Usage: miniTopSim.py tend dt")
    exit(0)

tend = float(sys.argv[1])
dt = float(sys.argv[2])

time = 0
    
surface = Surface()   

surface.write('basic_{}_{}.srf'.format(int(tend), int(dt)), time)    
     
while time < tend:
        
    adv.advance(surface,dt)
    dt = adv.timestep(dt, time, tend)
    time += dt
    surface.write('basic_{}_{}.srf'.format(int(tend), int(dt)), time)
    
surface.plot(tend, dt)