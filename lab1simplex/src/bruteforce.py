from src.utils import NpCanonicalForm, binomial_grid, subset_by_index
import numpy as np


def find_all_svs(cf: NpCanonicalForm):
    binom_table = binomial_grid(cf.n, cf.m)
    N = np.array(list(range(cf.n)))
    for idx in range(binom_table[-1, -1]):
        Nk = subset_by_index(N, binom_table, idx)
        AMNk = np.array([cf.A[:, i] for i in Nk]).T
        if np.linalg.det(AMNk) != 0:
            xNk = np.matmul(np.linalg.inv(AMNk), cf.b)
            xN = np.zeros(cf.n)
            for i in range(len(Nk)):
                xN[Nk[i]] = xNk[i]
            yield xN


def bruteforce(cf: NpCanonicalForm) -> np.array:
    """
    Алгоритм перебора решения задачи линейного программирования
    :param cf: параметры задачи в канонической форме
    :return: опорный вектор, минимизирующий целевую функцию c^T*x
    """
    svs = list(find_all_svs(cf))
    if not svs:  # опорные векторы не найдены, решения нет
        return np.zeros(cf.n)
    sv_min = svs[0]
    for sv in svs[1:]:
        if(np.dot(sv, cf.c) < np.dot(sv_min, cf.c)):
            sv_min = sv
    return sv_min
