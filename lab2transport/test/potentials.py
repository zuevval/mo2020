import numpy as np

import sys
sys.path.append('..')  # now environment can see package `src`
from src import potentials


def test_solve_transportation_potentials():
    # credit: Konykhovskiy P. V. (Конюховский П. В).
    a = np.array([27, 20, 43])
    b = np.array([33, 13, 27, 17])
    c_list = [(14, 28, 21, 28),
              (10, 17, 15, 24),
              (14, 30, 25, 21)]
    c = np.array([np.array(row) for row in c_list])
    x0_list = [(27,  None, None, None),
               (6,   13,   1,    None),
               (None, None, 26,  17)]
    x0 = np.array([np.array(row) for row in x0_list])

    u, v = potentials.calc_u_v(c, x0)
    v_expected = np.array([0, 7, 5, 1])
    u_expected = np.array([-14, -10, -20])
    assert (u_expected == u).all() and (v_expected == v).all()

    potentials.solve_transportation_potentials(a, b, c, x0)


if __name__ == "__main__":
    test_solve_transportation_potentials()
