

def test_raster_scan():
    import os
    import sys

    dirname = os.path.dirname(os.path.abspath(__file__))
    main_dir = os.path.join(dirname, os.path.pardir, os.path.pardir)
    code_dir = os.path.join(main_dir, "code")
    sys.path.insert(0, code_dir)

    import miniTopSim
    import surface as surf
    import parameters as par

    config_file = os.path.join(dirname, 'raster_scan_1.cfg')
    par.set_parameters(config_file)

    '''filename,'''
    sim = miniTopSim.simulate(config_file, True)
    surface_save = surf.load(os.path.join(dirname, 'raster_scan_1.srf_save'))

    assert sim.distance(surface_save) < 0.00292 