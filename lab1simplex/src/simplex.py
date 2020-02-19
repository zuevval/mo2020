import numpy as np
import scipy.special
from src.utils import NpCanonicalForm


def split_xkN(xkN:np.array) -> (np.array, np.array, np.array, np.array):
    Nk_plus = np.array(list(filter(lambda i: xkN[i] > 0, range(len(xkN)))))  # indices of positive components in xkN
    Nk0 = np.array(list(filter(lambda i: xkN[i] == 0, range(len(xkN))))) # indices of zero components in xkN
    xkNk_plus = np.array(list(filter(lambda x: x > 0, xkN)))
    xkNk0 = np.array(list(filter(lambda x: x == 0, xkN)))
    return xkNk0, xkNk_plus, Nk0.astype(int), Nk_plus.astype(int)


def new_AMNk(AMN: np.array, xkN:np.array, n_change:int) -> (np.array, np.array, np.array):
    # TODO implement for n_change > 0 (+ check determinant)
    m = len(AMN)
    xkNk0, xkNk_plus, Nk0, Nk_plus = split_xkN(xkN)
    Nk = Nk_plus
    Lk = Nk0
    Nk = np.sort(np.append(Lk[0:m - len(Nk_plus)], Nk))
    Lk = np.delete(Lk, range(m-len(Nk_plus)))
    AMNk = np.array([AMN[:, i] for i in Nk]).T
    return AMNk, Nk.astype(int), Lk.astype(int)


def calc_BNkM(AMNk):
    # TODO implement simplier method from book
    return np.linalg.inv(AMNk)


def simplex_step(cf: NpCanonicalForm, xkN: np.array) -> np.array:
    AMN, cN = cf.A, cf.c
    for AMNk_n_change in range(int(scipy.special.binom(cf.n, cf.m))):
        AMNk, Nk, Lk = new_AMNk(AMN, xkN, AMNk_n_change)
        BNkM = calc_BNkM(AMNk)
        cNk = np.array([cN[i] for i in Nk])
        ykM = np.matmul(BNkM.T, cNk)
        dkN = cN - np.matmul(AMN.T, ykM)
        dkLk = np.array([dkN[int(i)] for i in Lk])
        if np.min(dkLk) >= 0: # xkN is a solution
            print("solution found!")
            return xkN
        jk = list(filter(lambda j: dkLk[j] < 0, range(len(Lk))))[0]  # index of first negative component
        xkNk0, xkNk_plus, Nk0, Nk_plus = split_xkN(xkN)
        ukNk = np.matmul(BNkM, AMN[:,jk])
        xkNk = np.array([xkN[i] for i in Nk])
        if len(Nk_plus) == len(Nk) or max([ukNk[i] for i in filter(lambda j: j not in Nk_plus, Nk)]) < 0:
            ukN = [ukNk[list(Nk).index(i)] if i in Nk else 0 for i in range(cf.n)]
            ukN[jk] = -1
            theta_k = min([xkN[i]/ukN[i] for i in filter(lambda j: ukN[j] > 0, Nk)])
            return xkN - np.multiply(theta_k, ukN)
    return xkN
