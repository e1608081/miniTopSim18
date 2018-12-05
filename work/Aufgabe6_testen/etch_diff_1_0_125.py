import os
import sys

dirname = os.path.dirname(os.path.abspath(__file__))
main_dir = os.path.join(dirname, os.path.pardir, os.path.pardir)
code_dir = os.path.join(main_dir, "code")
sys.path.insert(0, code_dir)
import plot

plot.plot(os.path.join(dirname, 'etch_dx0_125_10_0.srf'),
		  os.path.join(dirname, 'etch_dx1_10_1.srf'))