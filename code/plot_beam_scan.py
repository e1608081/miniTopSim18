import sys
import time

import advance as adv
import parameters as par
import matplotlib.pyplot as plt
import miniTopSim
import plot
from sputtering import init_sputtering
from surface import Surface



def simulate(config_file, setconfig=True):
    """
    Performs a simulation with the parameters given inside the config file

        :param config_file: parameters for simulation

        :return the surface object after last simulation step
    """
    
    if setconfig is True:
        try:
            par.set_parameters(config_file)
        except FileNotFoundError as err:
            print(err)
            sys.exit(-1)


    config = "config"
    if (config_file.find(".cfg") != -1):
        config = config_file.replace(".cfg", "")


    tend = float(par.TOTAL_TIME)
    dt = float(par.TIME_STEP)
    
    time = 0

    init_sputtering()
    surface = Surface()
    surface_filename = '{}.srf'.format(config)
    surface.write(surface_filename, time, 'w')

    while time < tend:
#        print("time",time)
#        print("tend",tend)
#        print("dt",dt)
        #print("dwell",par.DWELL_TIME)
        #print("Beam dwell in main", beam.dwell_time)
       
        adv.advance(surface, dt)        
        dt = adv.timestep(dt, time, tend)
        time += dt
        surface.write(surface_filename, time, 'a')

    return surface_filename, surface

if __name__ == "__main__":

    time_start = time.process_time()
    print("This process can take some time depending on config")

    filename_scan, surface_scan = simulate(sys.argv[1])
    print("Surface done")
    time_used = time.process_time() - time_start
    print("Computing time = " + str(time_used))

    srf_scan = plot.loadFile(filename_scan)


    surface_scan = srf_scan[-1]
    xValues_scan = surface_scan["xVals"]
    yValues_scan = surface_scan["yVals"]
    plt.grid(True, 'major')
    plt.xlabel('x-values in nm')
    plt.ylabel('y-values in nm')


    plt.plot(xValues_scan, yValues_scan, label='Scanned Beam')
    plt.legend()
    plt.show()
