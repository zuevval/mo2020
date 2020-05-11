from typing import Callable

from cutting_plane import cutting_plane_alg as cut_alg
from convert_utils import LinearProblem
import numpy as np
import logging
import sys


def example():
    A = np.array([
        [0, 1],
        [0, -1],
        [1, 0],
        [-1, 0]
    ])
    b = np.array([1, 1, 1, 1])
    c = np.array([-1, -1])
    phi: Callable[[np.array], float] = lambda x: x[0] ** 2 + x[1] ** 2 - 1
    phi_subgrad: Callable[[np.array], np.array] \
        = lambda x: np.array([2 * x[0], 2 * x[1]])
    # minimum in (1/sqrt(2), 1/sqrt(2))
    lp = LinearProblem(np.array(A), np.array(b), np.array(c))
    res_x = cut_alg(lp, phi, phi_subgrad)
    logging.info(res_x)


if __name__ == "__main__":
    logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
    example()

