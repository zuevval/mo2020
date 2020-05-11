from typing import Tuple

import numpy as np  # type: ignore
import logging
from utils import binomial_grid, subset_by_index
from convert_utils import CanonicalForm


def split_xkN(xkN: np.array) -> Tuple[np.array, np.array, np.array, np.array]:
    Nk_plus = np.array(list(filter(lambda i: xkN[i] > 0, range(
        len(xkN)))))  # indices of positive components in xkN
    Nk0 = np.array(list(filter(lambda i: xkN[i] == 0, range(
        len(xkN)))))  # indices of zero components in xkN
    xkNk_plus = np.array(list(filter(lambda x: x > 0, xkN)))
    xkNk0 = np.array(list(filter(lambda x: x == 0, xkN)))
    return xkNk0, xkNk_plus, Nk0.astype(int), Nk_plus.astype(int)


def add_to_Nk(Nk0: np.array, binom_table: np.ndarray,
              binom_index: int) -> np.array:
    """
    Augments indices of positive basis vectors so that they form
    a square matrix when put together
    :param Nk0: an array of indices of positive vectors
    :param binom_table: binomial table created by `utils.binomial_grid` method
    :param binom_index: iteration of basis change (0 to C[|Nk_plus|,m-|Nk0|])
    :return: indices from `Nk0` that should be added to `Nk_plus` in order to
    make a square matrix from vectors with those indices
    """
    return subset_by_index(Nk0, binom_table, binom_index)


def new_AMNk(AMN: np.array, xkN: np.array, binom_table: np.ndarray,
             binom_idx: int) -> Tuple[np.array, np.array, np.array]:
    xkNk0, xkNk_plus, Nk0, Nk_plus = split_xkN(xkN)
    Nk = np.sort(
        np.append(Nk_plus, add_to_Nk(Nk0, binom_table, binom_idx))).astype(int)
    Lk = np.array(list(filter(lambda idx: idx not in Nk, Nk0))).astype(int)
    AMNk = np.array([AMN[:, i] for i in Nk]).T
    return AMNk, Nk.astype(int), Lk.astype(int)


class nextBParams:
    def __init__(self, BNkM: np.ndarray, ukNk: np.array, ik: int):
        self.m = len(ukNk)
        dim_correct = len(BNkM) == self.m
        for row in BNkM:
            dim_correct &= len(row) == self.m
        dim_correct &= ik < self.m
        if not dim_correct:
            raise Exception("dimensions mismatch")
        self.BNkM, self.ukNk, self.ik = BNkM, ukNk, ik


def calc_BNkM(AMNk: np.ndarray) -> np.ndarray:
    """
    Rejected calculation with usage of previous step info because it requires
    that sets of indices of appended vectors differ only in one index at once,
    while `subset_by_index(...)` works differently
    :param AMNk: a square matrix
    :return: an inverse matrix
    """
    return np.linalg.inv(AMNk)


def starting_vector(cf: CanonicalForm) -> np.array:
    binom_table = binomial_grid(cf.n, cf.m)
    N = np.array(list(range(cf.n)))
    for idx in range(binom_table[-1, -1]):
        Nk = subset_by_index(N, binom_table, idx)
        AMNk = np.array([cf.A[:, i] for i in Nk]).T
        if np.linalg.det(AMNk) != 0:
            xNk = np.matmul(np.linalg.inv(AMNk), cf.b)
            if np.min(xNk) < 0:
                continue  # we need only vectors with all positive components
            xN = np.zeros(cf.n)
            for i in range(len(Nk)):
                xN[Nk[i]] = xNk[i]
            return xN
    logging.info("error, initial vector not found!")
    return np.zeros(cf.n)


def simplex_step(cf: CanonicalForm, xkN: np.array) -> Tuple[
        np.array, np.array, bool]:
    """
    Performs one step of simplex algorithm.
    Parameters names corresponds to those in a book by Petuhov (Петухов) et al,
    p. 88. Algorithm is based on the same source.
    :param cf: parameters of problem in canonical form
    :param xkN: a starting vector: any support vector for a set defined by `cf`
    :return:
        1. The next approximation - also a support vector (the value of target
        function is not increased)
        2. Nk - a list of basis vectors of this support vector in A[M,N]
        3. A boolean equal to True if iterations are to be stopped,
         otherwise False
    """
    AMN, cN = cf.A, cf.c
    _, _, Nk0, Nk_plus = split_xkN(
        xkN)  # Nk0 - indices of zero components of xk, Nk+ - positive
    binomGrid = binomial_grid(len(Nk0), cf.m - len(
        Nk_plus))  # auxiliary structure for building combinations that are
    # added to A[M,Nk+]
    for binom_idx in range(binomGrid[-1, -1]):
        # augment Nk+ to Nk so that A[M, Nk] is square
        AMNk, Nk, Lk = new_AMNk(AMN, xkN, binomGrid, binom_idx)
        # if determinant of the square matrix is 0, continue
        if np.linalg.det(AMNk) == 0:
            continue
        BNkM = calc_BNkM(AMNk)  # calculating inverse matrix for A[M, Nk]
        cNk = np.array([cN[i] for i in Nk])
        ykM = np.matmul(BNkM.T, cNk)
        dkN = cN - np.matmul(AMN.T, ykM)
        dkLk = np.array([dkN[int(i)] for i in Lk])
        if np.min(dkLk) >= -1e-4:  # xkN is already the optimal vector
            logging.info("solution found at iteration " + str(binom_idx + 1))
            return xkN, Nk, True
        # index of first negative components in dkLk
        jk = Lk[
            list(filter(lambda j: dkLk[j] < -1e-4, range(len(Lk))))[0]]
        xkNk0, xkNk_plus, Nk0, Nk_plus = split_xkN(xkN)
        ukNk = np.matmul(BNkM, AMN[:, jk])
        if np.max(ukNk) <= -1e-4:  # target function is not lower bounded
            logging.info("solution does not exist")
            return np.array([np.inf for _ in range(cf.n)]), Nk, True
        if len(Nk_plus) == len(Nk) or max(
                [ukNk[i] for i in
                 filter(lambda j: Nk[j] not in Nk_plus, range(len(Nk)))]) < 0:
            ukN = [ukNk[list(Nk).index(i)] if i in Nk else 0 for i in
                   range(cf.n)]
            ukN[jk] = -1
            theta_k = min(
                [xkN[i] / ukN[i] for i in filter(lambda j: ukN[j] > 0, Nk)])
            return xkN - np.multiply(theta_k, ukN), Nk, False
    logging.info("iterations ended with no result, something went wrong.")
    return xkN, np.array([]), True  # something went wrong


class SimplexResult:
    def __init__(self, x: np.array, Nk: np.array, is_solution: bool):
        self.x, self.Nk, self.is_solution = x, Nk, is_solution


def simplex_alg(cf: CanonicalForm, x_start: np.array,
                   max_iter: int = 20) -> SimplexResult:
    xkN, Nk = x_start, np.array([])
    for _ in range(max_iter):
        xkN, Nk, stopIteration = simplex_step(cf, xkN)
        if stopIteration:
            break
    if len(Nk) == 0 or np.max(xkN) == np.inf:
        return SimplexResult(xkN, Nk, False)
    else:
        return SimplexResult(xkN, Nk, True)
