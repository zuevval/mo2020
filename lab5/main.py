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


def our_problem():
    def phi(x: np.array) -> float:
        phi1 = 3*x[0]**2 + x[1]**2 - 1
        phi2 = x[0]**2 + (x[1] - 0.5)**2 - 0.5
        phi3 = 3*x[0]**2 + x[1]**2-x[2]
        return np.array([phi1, phi2, phi3])

    def grad_phi(x: np.array, i: int) -> float:
        if i == 1:
            df_x = 6*x[0] # производная по x
            df_y = 2*x[1] # производная по y
        elif i == 2:
            df_x = 2*x[0] # производная по x
            df_y = 2*x[1] - 1 # производная по y
        elif i == 3:
            df_x = 2*x[0] # производная по x
            df_y = 2*x[1] # производная по y
            df_z = -1 # производная по z
            return sqrt(pow(df_x,2)+pow(df_y,2) + pow(df_z,2))

        return sqrt(pow(df_x,2)+pow(df_y,2))

    def I(x: np.array, m: int) -> list:
        I = []
        i = 0
        phi_values = phi(x)
        maxPhi = np.max(phi_values[0], phi_values[1], phi_values[2])
        while i < m:
            if phi_values[i] == maxPhi:
                I.append(i + 1) # для наглядности i + 1, чтоб номера индексов совпали с номером ограничения 
            i = i + 1
        return I

    def phi_subgrad(x: np.array) -> float:
        m = 3 # количество ограничений 
        I_array = I(x,m)
        return grad_phi(x, I_array[0]) # можно любой индекс, возьмем первый

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
    #example()
    logging.info("-------- our problem ---------")
    our_problem()

