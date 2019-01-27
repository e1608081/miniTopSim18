# -*- coding: utf-8 -*-
"""
Main file for miniTopSim
"""

import sys
import os
#from time import clock
from time import perf_counter

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
    
    #time_start = clock()
    time_start = perf_counter()
    
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
    surface_filename = '{}.srf'.format(config)
    surface.write(surface_filename, time, 'w')

    while time < tend:
        try:
            adv.advance(surface,dt)
        except ValueError as Err:
            print(Err)
            break
        else:
            dt = adv.timestep(dt, time, tend)
            time += dt
        finally:
            surface.write(surface_filename, time, 'a')
    
    #time_compute = clock() - time_start
    time_compute = perf_counter() - time_start

    print("Computing Time = " + str(time_compute))

    if par.PLOT_SURFACE and do_plotting:
        if os.path.isfile(surface_filename + "_save"):
            plot.plot(surface_filename, surface_filename + "_save")
        else:
            plot.plot(surface_filename)

    return surface


if __name__ == "__main__":
#    try:
        simulate(sys.argv[1])
#    except IndexError:
#        print('No file specified as systemargument')
#        sys.exit(-1)
