import os

import plot_beam
import matplotlib.pyplot as plt
import surface as surf
import parameters as par
import plot

config_file = 'erf.cfg'
par.set_parameters(config_file)

par.INITIAL_SURFACE_TYPE = 'Flat'
par.TOTAL_TIME = 1000
par.BEAM_CURRENT = 1e-12
par.XMAX = 100
par.XMIN = -100
filename, surface = plot_beam.simulate(config_file, False)

srf = plot.loadFile(filename)
if os.path.isfile(filename + '_save'):
    srf_save = plot.loadFile(filename + '_save')

else:
    raise FileNotFoundError(filename + '_save not found')

srf_save = plot.loadFile(filename + '_save')
surface_save = surf.load(filename + '_save')
distance = surface.distance(surface_save)
print('Distance between surface is ' + str(distance))
# assert distance < 0.00292933136297

surface_final = srf[-1]
surface_final_save = srf_save[-1]
xValues = surface_final["xVals"]
yValues = surface_final["yVals"]
xValues_save = surface_final_save["xVals"]
yValues_save = surface_final_save["yVals"]
plt.grid(True, 'major')
plt.xlabel('x-values in nm')
plt.ylabel('y-values in nm')

plt.plot(xValues, yValues, label='new calculation')
plt.plot(xValues_save, yValues_save, label='saved file')
plt.legend()
plt.show()
