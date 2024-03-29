"""
Test script that compares newly simulated surface with it's backup
"""


def test_cosine():
    import os
    import sys

    dirname = os.path.dirname(os.path.abspath(__file__))
    main_dir = os.path.join(dirname, os.path.pardir, os.path.pardir)
    code_dir = os.path.join(main_dir, "code")
    sys.path.insert(0, code_dir)
    import miniTopSim
    import surface
    import plot

    simulated = miniTopSim.simulate(os.path.join(dirname, 'cosine.cfg'), False)
    saved = surface.load(os.path.join(dirname, 'cosine.srf_save'))

    assert simulated.distance(saved) < 0.001
