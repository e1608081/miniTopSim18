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
    
    simulated = miniTopSim.simulate(os.path.join(code_dir, 'etch_dx1.cfg'), False)
    saved = surface.load(os.path.join(code_dir, 'etch_dx1_10_1.srf_save'))
    
    assert simulated.distance(saved) < 0.00215
