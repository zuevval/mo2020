import copy
import logging

from convert_utils import LinearProblem, CanonicalForm,\
    dual_problem, back_to_primal
from typing import Callable
import numpy as np
import bruteforce
from simplex import simplex_alg


class CuttingPlaneData:
    def __init__(self, xk: np.array, yk:np.array,
                 phi_subgrad: Callable[[np.array], np.array],
                 phi: Callable[[np.array], float], lp: LinearProblem):
        self.xk, self.yk, self.phi_subgrad, self.phi, self.lp \
            = xk, yk, phi_subgrad, phi, lp
        self.if_in_omega: Callable[[np.array], bool] \
            = lambda x: phi(x) <= 0


def starting_stage(lp: LinearProblem) -> [np.array, np.array]:
    cf: CanonicalForm = dual_problem(lp)
    lin_solution = bruteforce.bruteforce(cf)
    y0, Nk = lin_solution.x, lin_solution.Nk
    x0 = back_to_primal(Nk=Nk, primal=lp)
    return x0, y0




def cutting_plane_iteration(data: CuttingPlaneData)\
        -> CuttingPlaneData:
    # find a_[k+1], b_[k+1]
    a_k1: np.array = data.phi_subgrad(data.xk)
    b_k1: float = -data.phi(data.xk) + np.dot(a_k1, data.xk)

    # S_k -> S_[k+1]: append row to matrix and number to vector
    A_next = np.concatenate([data.lp.A, np.array([a_k1])])
    b_next = np.append(data.lp.b, b_k1)
    logging.debug("added plane: " + str(a_k1) + " * x  <= " + str(b_k1))
    lp_next = LinearProblem(A=A_next, b=b_next, c=data.lp.c)

    # solve linear programming problem
    cf_dual: CanonicalForm = dual_problem(lp_next)
    #lin_res = simplex_alg(cf_dual, data.yk)
    #y_k1, Nk = lin_res.x, lin_res.Nk
    #inv_AMNk = np.linalg.inv(cf_dual.A[:, Nk])
    lin_res = bruteforce.bruteforce(cf_dual) #  alternative to lines above
    assert (abs(cf_dual.A.dot(lin_res.x) - cf_dual.b) < 1e-3).all()
    y_k1, Nk = lin_res.x, lin_res.Nk
    x_k1 = back_to_primal(Nk=Nk, primal=lp_next)
    logging.debug(A_next[Nk].dot(x_k1) - b_next[Nk])
    # assert (A_next.dot(x_k1) <= b_next).all() # TODO
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
    assert (lp.A.dot(x0) <= lp.b).all()
    data = CuttingPlaneData(xk=x0, yk=y0, phi_subgrad=phi_subgrad,
                            phi=phi, lp=lp)
    for _ in range(max_iter):
        logging.debug("xk: " + str(data.xk))
        # step 1 (see presentation in Teams, lecture 5, slide 5)
        if data.if_in_omega(data.xk):
            logging.debug("x_k in omega")
            np.savetxt("r-vis/data/A.txt", data.lp.A)
            np.savetxt("r-vis/data/b.txt", data.lp.b)
            return data.xk
        # step 2
        data_next = cutting_plane_iteration(copy.deepcopy(data))
        if np.linalg.norm(data_next.xk - data.xk) < eps:
            logging.debug("x_{k+1} and x_k are close enough to stop")
            np.savetxt("r-vis/data/A.txt", data_next.lp.A)
            np.savetxt("r-vis/data/b.txt", data_next.lp.b)
            return data_next.xk
        data = data_next
    # the code below is normally unreachable
    logging.warning("cutting-plane method "
                    "reached maximal number of iterations")
    return data.xk