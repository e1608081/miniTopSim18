# -*- coding: utf-8 -*-
"""
Main file for miniTopSim
"""

import sys

import advance as adv
import parameters as par
import plot
from sputtering import init_sputtering
from surface import Surface

try:
    par.set_parameters(sys.argv[1])
except FileNotFoundError as err:
    print(err)
    sys.exit(-1)
except IndexError:
    print('No file specified as systemargument')
    sys.exit(-1)

config = "config"
if sys.argv[1].find(".cfg") != -1:
    config = sys.argv[1].replace(".cfg", "")

tend = float(par.TOTAL_TIME)
dt = float(par.TIME_STEP)

time = 0

init_sputtering()
surface = Surface()

surface.write('{}_{}_{}.srf'.format(config, int(tend), int(dt)), time, 'w')

while time < tend:

    adv.advance(surface, dt)
    dt = adv.timestep(dt, time, tend)
    time += dt
    surface.write('{}_{}_{}.srf'.format(config, int(tend), int(dt)), time, 'a')

if par.PLOT_SURFACE:
    plot.plot('{}_{}_{}.srf'.format(config, int(tend), int(dt)))

