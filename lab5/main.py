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


def our_problem():
    def phi(x: np.array) -> float:
        phi1 = 3*x[0]**2 + x[1]**2 - 1
        phi2 = x[0]**2 + (x[1] - 0.5)**2 - 0.5
        phi3 = 3*x[0]**2 + x[1]**2-x[2]
        return np.max([phi1, phi2, phi3])

    def phi_subgrad(x: np.array) -> np.array:
        pass  # TODO

    A = np.array([
        [1, 0, 0],
        [-1, 0, 0],
        [0, 1, 0],
        [0, -1, 0],
        [0, 0, 1],
        [0, 0, -1]
    ])
    b = np.array([1/np.sqrt(3), 1/np.sqrt(3), 1,
                  (1 - np.sqrt(2))/2, -1, 0])
    c = np.array([0, 0, 1])
    lp = LinearProblem(A, b, c)
    res_x = cut_alg(lp, phi, phi_subgrad)
    logging.info(res_x)


if __name__ == "__main__":
    logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
    logging.info("------- simple example -------")
    example()
    logging.info("-------- our problem ---------")
    our_problem()

