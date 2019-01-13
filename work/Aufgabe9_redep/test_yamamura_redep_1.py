"""
Test script that compares newly simulated surface with it's backup
"""


def test_yamamura_redep_1():
    import os
    import sys

    dirname = os.path.dirname(os.path.abspath(__file__))
    main_dir = os.path.join(dirname, os.path.pardir, os.path.pardir)
    code_dir = os.path.join(main_dir, "code")
    sys.path.insert(0, code_dir)
    import miniTopSim
    import surface
    import plot

    simulated = miniTopSim.simulate(os.path.join(dirname, 'yamamura_redep_1.cfg'), False)
    saved = surface.load(os.path.join(dirname, 'yamamura_redep.srf_save2'))

    assert simulated.distance(saved) < 0.0027584
