import numpy as np

def normal(x, y):
    """Calc normal vector and return it."""
    x = np.concatenate(([x[0]], x, [x[-1]]))
    y = np.concatenate(([y[0]], y, [y[-1]]))
    dx = _calc_vector(x)
    dy = _calc_vector(y)
    length = np.linalg.norm([dx, dy], axis=0)
    return dy / length, -dx / length

def _calc_vector(value):
    """Subtract end coordinate with start coordinate.

    :param value: coordinate array
    """
    delta = value[2:] - value[:-2]
    return delta

def calc_viewfactor(x, y):
    """Calculates the view-factor from surface parameters

    :return: nxn matrix representing the view-factor
    """
    # create arrays for nodes i & j (j has equal rows, i has equal columns)
    x_j = np.ones(shape=(x.size, x.size)) * x
    y_j = np.ones(shape=(y.size, y.size)) * y

    x_i = x_j.transpose()
    y_i = y_j.transpose()

    # calculate distances between nodes i & j
    x_ij = x_i - x_j
    y_ij = y_i - y_j

    # create arrays for the normal vectors of nodes i & j (j is transposed)
    x_normal, y_normal = normal(x, y)
    x_i_normal = np.ones_like(x_i) * x_normal
    y_i_normal = np.ones_like(y_i) * y_normal
    x_j_normal = np.ones_like(x_j) * x_normal[:, np.newaxis]
    y_j_normal = np.ones_like(y_j) * y_normal[:, np.newaxis]

    # calculate cosines of angles (cos_a, cos_b) with [(a*b)/(|a|*|b|)]
    # deactivate warnings when division by zero -> is covered with np.nan_to_num
    with np.errstate(divide='ignore', invalid='ignore'):
        cos_a = np.nan_to_num((x_j_normal * x_ij + y_j_normal * y_ij) / (
        np.sqrt(x_j_normal ** 2 + y_j_normal ** 2) * np.sqrt(x_ij ** 2 + y_ij ** 2)))
        cos_b = np.nan_to_num((x_i_normal * x_ij + y_i_normal * y_ij) / (
        np.sqrt(x_i_normal ** 2 + y_i_normal ** 2) * np.sqrt(x_ij ** 2 + y_ij ** 2)))

    print(cos_a)
    print(cos_b)

    # calculate distances between nodes i & j (d_ij)
    d_ij = np.sqrt(x_ij ** 2 + y_ij ** 2)

    # calculate surface length of node (delta_l)
    x_ext = np.concatenate(([2 * x[0] - x[1]], x, [3 * x[-1] - x[-2]]))
    y_ext = np.concatenate(([y[0]], y, [y[-1]]))
    d_x = x_ext[1:] - x_ext[:-1]
    d_y = y_ext[1:] - y_ext[:-1]
    d_l = np.sqrt(d_x ** 2 + d_y ** 2)
    d_l_avg = (d_l[1:] + d_l[:-1]) / 2
    delta_l = np.ones_like(x_i) * d_l_avg

    # calculate view-factor and mask out all values where cos_a<0 and cos_b<0 and all elements on diagonal
    # deactivate warnings when division by zero ->is covered with np.nan_to_num
    with np.errstate(divide='ignore', invalid='ignore'):
        f_ij = np.nan_to_num((cos_a * cos_b * delta_l) / (2 * d_ij))

    mask = (cos_a > np.zeros_like(x_i)) * (cos_b > np.zeros_like(x_i)) * np.invert(np.eye(x.size, dtype=bool))
    print(f_ij)
    print(mask)
    return f_ij*mask

x = np.arange(5)
y = np.array([0,0,-2,0,0])

print(x)
print(y)
print(calc_viewfactor(x,y))
