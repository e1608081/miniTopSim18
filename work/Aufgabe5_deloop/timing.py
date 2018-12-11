# -*- coding: utf-8 -*-
"""
plotscript for the computingtime
based on miniTopSim.py
"""


import sys
from time import clock
import matplotlib.pyplot as plt
import numpy as np

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

config = "config"
if(sys.argv[1].find(".cfg") != -1):
    config = sys.argv[1].replace(".cfg","")
    
tend = float(par.TOTAL_TIME)
dt = float(par.TIME_STEP)

#arrays for the plot
times = []
n = []


for i in np.fromfunction(lambda i: 10**(i/10-1), (11,)) :  
    
    dx = i  
    
    time_start = clock()
    
    time = 0
        
    surface = Surface()
    surface.set_dx(dx)   
    
    #surface.write('{}_{}_{}.srf'.format(config, int(tend), int(dt)), time, 'w')    
         
    while time < tend:
            
        adv.advance(surface,dt)
        dt = adv.timestep(dt, time, tend)
        time += dt
        #surface.write('{}_{}_{}.srf'.format(config, int(tend), int(dt)), time, 'a')
        
    time_program = clock() - time_start
    times.append(time_program)
    n.append(100/dx)
    print("x values = " + str(100/dx))
    print("Computingn Time = " + str(time_program))

    if par.PLOT_SURFACE:
        plot.plot('{}_{}_{}.srf'.format(config, int(tend), int(dt)))

plt.title("Computing time")
plt.plot(n, times, 'o-')
plt.ylabel("time [s]")
plt.xlabel("x-values [1]")
plt.yscale('log')
plt.xscale('log')
plt.grid()
plt.show()
        
