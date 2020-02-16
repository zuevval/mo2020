import numpy as np
from src.utils import NpCanonicalForm


def split_xkN(xkN:np.array) -> (np.array, np.array, np.array, np.array):
    Nk_plus = np.array(list(filter(lambda i: xkN[i] > 0, range(len(xkN)))))  # indices of positive components in xkN
    Nk0 = np.array(list(filter(lambda i: xkN[i] == 0, range(len(xkN))))) # indices of zero components in xkN
    xkNk_plus = np.array(list(filter(lambda x: x > 0, xkN)))
    xkNk0 = np.array(list(filter(lambda x: x == 0, xkN)))
    return xkNk0, xkNk_plus, Nk0, Nk_plus


def new_AMNk(AMN: np.array, xkN:np.array, n_change:int) -> (np.array, np.array, np.array):
    # TODO implement for n_change > 0 (+ check determinant)
    m = len(AMN)
    xkNk0, xkNk_plus, Nk0, Nk_plus = split_xkN(xkN)
    Nk = Nk_plus
    Lk = Nk0
    Nk = np.sort(np.append(Lk[0:m - len(Nk_plus)], Nk))
    Lk = np.delete(Lk, range(m-len(Nk_plus)))
    AMNk = np.array([AMN[:, int(i)] for i in Nk]).T
    return AMNk, Nk, Lk


def calc_BNkM(AMNk):
    # TODO implement simplier method from book
    return np.linalg.inv(AMNk)


def simplex_step(cf: NpCanonicalForm, xNk: np.array) -> np.array:
    # TODO implement
    return xNk
