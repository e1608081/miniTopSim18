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
        adv.advance(surface, dt)
        dt = adv.timestep(dt, time, tend)
        time += dt
        surface.write(surface_filename, time, 'a')

    return surface_filename, surface

if __name__ == "__main__":

    time_start = time.process_time()
    print("This process can take some time depending on config")
    filename_const, surface_const = simulate('constant.cfg')
    print("Constant surface done")
    filename_gauss, surface_gauss = simulate('gauss.cfg')
    print("Gaussian surface done")
    filename_erf, surface_erf = simulate('erf.cfg')
    print("Erf surface done")
    time_used = time.process_time() - time_start
    print("Computing time = " + str(time_used))

    srf_const = plot.loadFile(filename_const)
    srf_gauss = plot.loadFile(filename_gauss)
    srf_erf = plot.loadFile(filename_erf)


    surface_gauss = srf_gauss[-1]
    xValues_gauss = surface_gauss["xVals"]
    yValues_gauss = surface_gauss["yVals"]
    surface_erf = srf_erf[-1]
    xValues_erf = surface_erf["xVals"]
    yValues_erf = surface_erf["yVals"]
    surface_const = srf_const[-1]
    xValues_const = surface_const["xVals"]
    yValues_const = surface_const["yVals"]
    plt.grid(True, 'major')
    plt.xlabel('x-values in nm')
    plt.ylabel('y-values in nm')

    plt.plot(xValues_gauss, yValues_gauss, label='Gaussian Beam')
    plt.plot(xValues_erf, yValues_erf, label='Erf Beam')
    plt.plot(xValues_const, yValues_const, label='Constant Beam')
    plt.legend()
    plt.show()
