# -*- coding: utf-8 -*-
"""
Main file for miniTopSim
"""

import sys
import os
from time import clock

import advance as adv
import parameters as par
import plot
from sputtering import init_sputtering
from surface import Surface

def simulate(config_file, do_plotting=True):
    """
    Performs a simulation with the parameters given inside the config file
    
        :param config_file: parameters for simulation
        
        :return the surface object after last simulation step
    """
    
    time_start = clock()
    
    try:
        par.set_parameters(config_file)
    except FileNotFoundError as err:
        print(err)
        sys.exit(-1)
    
    config = "config"
    if(config_file.find(".cfg") != -1):
        config = config_file.replace(".cfg","")
        
    tend = float(par.TOTAL_TIME)
    dt = float(par.TIME_STEP)
    time = 0
    
    init_sputtering()
    surface = Surface()
    surface_filename = '{}_{}_{}.srf'.format(config, int(tend), int(dt))
    surface.write(surface_filename, time, 'w')    
         
    while time < tend:
        adv.advance(surface,dt)
        dt = adv.timestep(dt, time, tend)
        time += dt
        surface.write(surface_filename, time, 'a')
    
    time_compute = clock() - time_start
    print("Computing Time = " + str(time_compute))
    
    if par.PLOT_SURFACE and do_plotting:
        if os.path.isfile(surface_filename + "_save"):
            plot.plot(surface_filename, surface_filename + "_save")
        else:
            plot.plot(surface_filename)
    
    return surface

if __name__ == "__main__":
    try:
        simulate(sys.argv[1])
    except IndexError:
        print('No file specified as systemargument')
        sys.exit(-1)
