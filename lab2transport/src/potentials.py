import numpy as np


def calc_u_v(c:np.ndarray, x0:np.ndarray) -> (np.array, np.array):
    (m, n) = c.shape
    u, v = np.zeros(n), np.zeros(m)
    # TODO calc u, v
    return u, v


def solve_transportation_potentials(a:np.array, b:np.array, c:np.ndarray, x0:np.ndarray):
    pass
