

def test_erf():
    import os
    import sys

    dirname = os.path.dirname(os.path.abspath(__file__))
    main_dir = os.path.join(dirname, os.path.pardir, os.path.pardir)
    code_dir = os.path.join(main_dir, "code")
    sys.path.insert(0, code_dir)

    import plot_beam
    import surface as surf
    import parameters as par

    config_file = os.path.join(dirname, 'erf.cfg')
    par.set_parameters(config_file)

    # cfg overwrites
    par.INITIAL_SURFACE_TYPE = 'Flat'
    par.TOTAL_TIME = 1
    # par.BEAM_CURRENT = 1e-12
    par.XMAX = 500
    par.XMIN = -500
    filename, surface = plot_beam.simulate(config_file, False)

    if os.path.isfile(filename + '_save'):
        surface_save = surf.load(filename + '_save')

    else:
        raise FileNotFoundError(filename + '_save not found')

    distance = surface.distance(surface_save)
    assert distance < 0.00292933136297