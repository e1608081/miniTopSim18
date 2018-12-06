# -*- coding: utf-8 -*-
"""
Main file for miniTopSim
"""


import sys

from surface import Surface
import advance  as adv

import parameters as par
import plot

def simulate(config_file):
    """
    Performs a simulation with the parameters given inside the config file
    
        :param config_file: parameters for simulation
        
        :return the surface object after last simulation step
    """
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
    
    surface = Surface()   
    surface.write('{}_{}_{}.srf'.format(config, int(tend), int(dt)), time, 'w')    
         
    while time < tend:
        adv.advance(surface,dt)
        dt = adv.timestep(dt, time, tend)
        time += dt
        surface.write('{}_{}_{}.srf'.format(config, int(tend), int(dt)), time, 'a')
    
    if par.PLOT_SURFACE:
        plot.plot('{}_{}_{}.srf'.format(config, int(tend), int(dt)))
    
    return surface

if __name__ == "__main__":
    try:
        simulate(sys.argv[1])
    except IndexError:
        print('No file specified as systemargument')
        sys.exit(-1)
