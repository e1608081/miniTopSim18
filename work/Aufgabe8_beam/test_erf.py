


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
    import plot

    config_file = 'erf.cfg'
    par.set_parameters(config_file)

    par.INITIAL_SURFACE_TYPE = 'Flat'
    # par.TOTAL_TIME = 1000
    # par.BEAM_CURRENT = 1e-12
    # par.XMAX = 100
    # par.XMIN = -100
    filename, surface = plot_beam.simulate(config_file, False)

    srf = plot.loadFile(filename)
    if os.path.isfile(filename + '_save'):
        srf_save = plot.loadFile(filename + '_save')

    else:
        raise FileNotFoundError(filename + '_save not found')

    surface_save = surf.load(filename + '_save')
    distance = surface.distance(surface_save)
    assert distance < 0.00292933136297
