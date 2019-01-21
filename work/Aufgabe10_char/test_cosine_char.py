"""
Test the characteristics time integration method by comparing to a saved surface
"""


def test_cosine_char():
    import os
    import sys
    filedir = os.path.dirname(os.path.abspath(__file__))
    codedir = os.path.join(filedir, '..', '..', 'code')
    sys.path.insert(0, codedir)

    import miniTopSim
    import surface

    saved = surface.load(os.path.join(filedir, 'step_char.srf_save'))
    simulated = miniTopSim.simulate(os.path.join(filedir, 'step_char.cfg'), False)

    assert simulated.distance(saved) < 0.00292933136297
