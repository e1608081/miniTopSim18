"""
Test the vertical time integration method and the initialization by
surface file. First, calculate only half the time steps, then load the
resulting .srf and calculate the second half of the time steps.
Finally compare the surface to a saved .srf_save file.
"""


def test_cosine_vert_2():
    import os
    import sys
    filedir = os.path.dirname(os.path.abspath(__file__))
    codedir = os.path.join(filedir, '..', '..', 'code')
    sys.path.insert(0, codedir)

    import miniTopSim
    import surface

    miniTopSim.simulate(os.path.join(filedir, 'cosine_vert_1.cfg'), False)
    simulated = miniTopSim.simulate(os.path.join(filedir, 'cosine_vert_2.cfg'), False)
    saved = surface.load(os.path.join(filedir, 'cosine_vert.srf_save'))

    assert simulated.distance(saved) < 0.00292933136297
