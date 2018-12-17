"""
Test script that compares newly simulated surface with it's backup
"""


def test_table_10_1():
    import os
    import sys

    dirname = os.path.dirname(os.path.abspath(__file__))
    main_dir = os.path.join(dirname, os.path.pardir, os.path.pardir)
    code_dir = os.path.join(main_dir, "code")
    sys.path.insert(0, code_dir)
    import miniTopSim
    import surface
    import plot

    simulated = miniTopSim.simulate(os.path.join(dirname, 'table.cfg'),
                                    False)
    saved = surface.load(os.path.join(dirname, 'table.srf_save'))

    assert simulated.distance(saved) < 0.0027584
