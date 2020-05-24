from typing import Callable

from cutting_plane import cutting_plane_alg as cut_alg
from convert_utils import LinearProblem
from math import sqrt
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


def pretty_print(our_res: np.array, phi: callable):
    our_res[2] += 0.7
    logging.info(our_res)
    logging.info("phi(x*): " + str(phi(our_res)))


def our_problem():
    def phi_components(x: np.array) -> np.array:
        phi1 = 3*x[0]**2 + x[1]**2 - 1
        phi2 = x[0]**2 + (x[1] - 0.5)**2 - 0.5
        phi3 = 3*x[0]**2 + x[1]**2 - x[2] - 1
        return np.array([phi1, phi2, phi3])

    def phi(x: np.array) -> float:
        return np.max(phi_components(x))

    def grad_phi(x: np.array, i: int) -> np.array:
        if i == 1:
            df_x = 6*x[0] # производная по x
            df_y = 2*x[1] # производная по y
        elif i == 2:
            df_x = 2*x[0] # производная по x
            df_y = 2*x[1] - 1 # производная по y
        elif i == 3:
            df_x = 6*x[0] # производная по x
            df_y = 2*x[1] # производная по y
            df_z = -1  # производная по z
            return np.array([df_x, df_y, df_z])

        return np.array([df_x, df_y, 0])

    def phi_subgrad(x: np.array) -> np.array:
        idx = np.argmax([phi_components(x)]) + 1
        return grad_phi(x, idx)

    A = np.array([
        [1, 0, 0],
        [-1, 0, 0],
        [0, 1, 0],
        [0, -1, 0],
        [0, 0, -1],
        [0, 0, 1]
    ])
    b = np.array([0.75, 0.75, 1,
                 (np.sqrt(2) - 1)/2 - 0.75, 1, 0])
    c = np.array([0, 0, 1])
    lp = LinearProblem(A, b, c)
    res_x = cut_alg(lp, phi, phi_subgrad, max_iter=15)
    pretty_print(res_x, phi)


if __name__ == "__main__":
    logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
    logging.info("\n\n------- simple example -------")
    example()
    logging.info("\n\n-------- our problem ---------")
    our_problem()

