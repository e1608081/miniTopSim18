import numpy as np
MAX_SEGLEN = 2.

x = np.array([2.1, 2.5, 4., 5., 5.5, 7., 13.])
y = np.array([0., 0., 0., 0., 0., 0., 0.])

x_dist2 = np.abs(x[2:] - x[:-2])
y_dist2 = np.abs(y[2:] - y[:-2])
dist2 = np.insert(np.sqrt(x_dist2 ** 2 + y_dist2 ** 2), (0, 1), 0)
remove_mask = np.where(dist2 < MAX_SEGLEN)