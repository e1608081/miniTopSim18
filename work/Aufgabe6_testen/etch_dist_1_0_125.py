import os
import sys

dirname = os.path.dirname(os.path.abspath(__file__))
main_dir = os.path.join(dirname, os.path.pardir, os.path.pardir)
code_dir = os.path.join(main_dir, "code")
sys.path.insert(0, code_dir)
import surface

surface1 = surface.load('etch_dx0_125.srf')
surface2 = surface.load('etch_dx1.srf')
print('Distance = ', surface1.distance(surface2))

