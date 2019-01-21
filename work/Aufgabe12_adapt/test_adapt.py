import os, sys

dirname = os.path.dirname(os.path.abspath(__file__))
filedir = os.path.dirname(__file__)
codedir = os.path.join(filedir, '..', '..', 'code')
sys.path.insert(0, codedir)


def test_adapt():
    import numpy as np
    import miniTopSim
    import surface

    srf = surface.load(os.path.join(dirname, 'temp.srf'))

    x = np.array(srf.x)
    y = np.array(srf.y)

    x_dist = x[1:]-x[:-1]
    y_dist = y[1:]-y[:-1]
    dist = np.sqrt(np.square(x_dist)+np.square(y_dist))
    assert(np.all(dist < 2.))

    normals_x, normals_y = np.array(srf.normal())
    angles = np.angle(normals_x + 1j*normals_y, deg=True)
    angles_diff = angles[1:]-angles[:-1]
    assert(np.all(angles_diff <= 10.))
