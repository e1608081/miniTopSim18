import numpy as np
import matplotlib.pyplot as plt

import surface

MAX_ANGLE = 10

Surface = surface.load('temp.srf')
vect_x, vect_y = Surface.normal()
angs = np.angle((vect_x+vect_y*1j), deg=True)
# insert 0 to match angle diff with point index
ang_diff = np.insert(np.abs(angs[1:]-angs[:-1]), 0, 0)
add_index = np.where(ang_diff > MAX_ANGLE)
plt.close()
plt.plot(Surface.x, Surface.y, '.-', Surface.x+vect_x, Surface.y+vect_y, 'r.')
plt.ylim((-11,2))
plt.xlim((-2,11))