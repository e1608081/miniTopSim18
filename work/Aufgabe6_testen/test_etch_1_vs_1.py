"""
Test Skript that compares newly simulated surface with it's backup
Should work
"""

def test_etch_1_vs_1():
    import os
    import sys
    
    dirname = os.path.dirname(os.path.abspath(__file__))
    main_dir = os.path.join(dirname, os.path.pardir, os.path.pardir)
    code_dir = os.path.join(main_dir, "code")
    sys.path.insert(0, code_dir)
    import miniTopSim
    import surface
    import plot
    
    simulated = miniTopSim.simulate(os.path.join(dirname, 'etch_dx1.cfg'), False)
    saved = surface.load(os.path.join(dirname, 'etch_dx1.srf_save'))
    
    #measured value: (formerly: 0.0135894565864)
    assert simulated.distance(saved) < 0.00292933136297 #(formerly: 0.0027584)
