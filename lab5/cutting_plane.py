import copy
import logging

from convert_utils import LinearProblem, CanonicalForm,\
    dual_problem, back_to_primal
from typing import Callable
import numpy as np
import bruteforce


class CuttingPlaneData:
    def __init__(self, xk: np.array, yk:np.array,
                 phi_subgrad: Callable[[np.array], np.array],
                 phi: Callable[[np.array], float], lp: LinearProblem):
        self.xk, self.phi_subgrad, self.phi, self.lp \
            = xk, phi_subgrad, self.phi, lp
        self.if_in_omega: Callable[[np.array], bool] \
            = lambda x: phi(x) <= 0


def starting_stage(lp: LinearProblem) -> [np.array, np.array]:
    cf: CanonicalForm = dual_problem(lp)
    y0 = bruteforce.bruteforce(cf)  # TODO solve with simplex method (?)
    x0 = back_to_primal(y0, lp, cf)  # TODO implement back_to_primal
    return x0, y0


def cutting_plane_iteration(data: CuttingPlaneData)\
        -> CuttingPlaneData:
    # find a_[k+1], b_[k+1]
    a_k1: np.array = data.phi_subgrad(data.xk)
    b_k1: float = -data.phi(data.xk) + np.dot(a_k1, data.xk)

    # S_k -> S_[k+1]: append row to matrix and number to vector
    A_next = np.concatenate([data.lp.A, np.array([a_k1])])
    b_next = np.append(data.lp.b, b_k1)
    lp_next = LinearProblem(A=A_next, b=b_next, c=data.lp.c)

    # solve linear programming problem
    cf_dual: CanonicalForm = dual_problem(lp_next)
    # TODO apply simplex method using yk instead of brute force
    y_k1 = bruteforce.bruteforce(cf_dual)
    x_k1 = back_to_primal(x=y_k1, primal=lp_next, dual=cf_dual)
    return CuttingPlaneData(xk=x_k1, yk=y_k1,
                            phi_subgrad=data.phi_subgrad,
                            phi=data.phi, lp=lp_next)


def cutting_plane_alg(
        lp: LinearProblem,
        phi: Callable[[np.array], float],
        phi_subgrad: Callable[[np.array], np.array],
        eps: float = 0.01,
        max_iter: int = 50) -> np.array:
    """
    Cutting plane minimization algorithm
    :param lp: linear problem with params A, b, c: cTx->min, x in S:={x|Ax<=b}
    :param phi: limitation: we search x in Omega := {x|phi(x)<=0},
    phi is a convex function and Omega: subset of S
    :param phi_subgrad: sub-gradient of phi,
    i. e. phi(x) >= phi(x0) + phi_subgrad(x0)(x-x0) for all x, x0 in S
    :param eps: stopping criteria
    :param max_iter: max number of iterations
    :return: solution of problem cTx->min, x in Omega
    """
    x0, y0 = starting_stage(lp)
    data = CuttingPlaneData(xk=x0, yk=y0, phi_subgrad=phi_subgrad,
                            phi=phi, lp=lp)
    for _ in range(max_iter):
        # step 1 (see presentation in Teams, lecture 5, slide 5)
        if data.if_in_omega(data.xk):
            return data.xk
        data_next = cutting_plane_iteration(copy.deepcopy(data))
        if np.norm(data_next.xk - data.xk) < eps:
            return data_next.xk
    logging.warning("cutting-plane method "
                    "reached maximal number of iterations")
    return np.array([np.inf for _ in range(len(data.xk))])