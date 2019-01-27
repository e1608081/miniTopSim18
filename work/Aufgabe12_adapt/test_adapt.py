def test_adapt():
    import os
    import sys

    dirname = os.path.dirname(os.path.abspath(__file__))
    main_dir = os.path.join(dirname, os.path.pardir, os.path.pardir)
    code_dir = os.path.join(main_dir, "code")
    sys.path.insert(0, code_dir)
    import numpy as np
    import surface

    srf = surface.load(os.path.join(dirname, 'cosine.srf_save'))

    x = np.array(srf.x)
    y = np.array(srf.y)

    x_dist = x[1:]-x[:-1]
    y_dist = y[1:]-y[:-1]
    dist = np.abs(np.sqrt(x_dist**2 + y_dist**2))
    assert np.all(dist <= 2.01)

    normals_x, normals_y = np.array(srf.normal())
    angles = np.angle(normals_x + 1j*normals_y, deg=True)
    angles_diff = angles[1:]-angles[:-1]
    assert np.all(angles_diff <= 10.)
