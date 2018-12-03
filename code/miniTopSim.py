# -*- coding: utf-8 -*-
"""
Main file for miniTopSim
"""


import sys

from surface import Surface
import advance  as adv

import parameters as par
import plot


try:
    par.set_parameters(sys.argv[1])
except FileNotFoundError as err:
    print(err)
    sys.exit(-1)
except IndexError:
    print('No file specified as systemargument')
    sys.exit(-1)

tend = float(par.TOTAL_TIME)
dt = float(par.TIME_STEP)

time = 0
    
surface = Surface()   

surface.write('basic_{}_{}.srf'.format(int(tend), int(dt)), time, 'w')    
     
while time < tend:
        
    adv.advance(surface,dt)
    dt = adv.timestep(dt, time, tend)
    time += dt
    surface.write('basic_{}_{}.srf'.format(int(tend), int(dt)), time, 'a')
  
if par.PLOT_SURFACE:
    plot.plot('basic_{}_{}.srf'.format(int(tend), int(dt)))