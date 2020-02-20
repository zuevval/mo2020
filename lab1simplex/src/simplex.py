import numpy as np
from src.utils import NpCanonicalForm, binomial_grid, subset_by_index


def split_xkN(xkN:np.array) -> (np.array, np.array, np.array, np.array):
    Nk_plus = np.array(list(filter(lambda i: xkN[i] > 0, range(len(xkN)))))  # indices of positive components in xkN
    Nk0 = np.array(list(filter(lambda i: xkN[i] == 0, range(len(xkN))))) # indices of zero components in xkN
    xkNk_plus = np.array(list(filter(lambda x: x > 0, xkN)))
    xkNk0 = np.array(list(filter(lambda x: x == 0, xkN)))
    return xkNk0, xkNk_plus, Nk0.astype(int), Nk_plus.astype(int)


def add_to_Nk(Nk0: np.array,  binom_table: np.ndarray, binom_index: int)->np.array:
    """
    Дополняет индексы положительных базисных векторов так, чтобы составленная из них матрица была квадратной
    :param Nk0: массив индексов положительных переменных
    :param binom_table: биномиальная таблица, созданная функцией `utils.binomial_grid`
    :param binom_index: номер итерации смены базиса от 0 до C[|Nk_plus|,m-|Nk0|]
    :return: индексы из `Nk0`, которые надо добавить к `Nk_plus`, чтобы получить квадратную матрицу
    """
    return subset_by_index(Nk0, binom_table, binom_index)


def new_AMNk(AMN: np.array, xkN:np.array, binom_table:np.ndarray, binom_idx:int) -> (np.array, np.array, np.array):
    m = len(AMN)
    xkNk0, xkNk_plus, Nk0, Nk_plus = split_xkN(xkN)
    Nk = np.sort(np.append(Nk_plus, add_to_Nk(Nk0, binom_table, binom_idx))).astype(int)
    Lk = np.array(list(filter(lambda idx: idx not in Nk, Nk0))).astype(int)
    AMNk = np.array([AMN[:, i] for i in Nk]).T
    return AMNk, Nk.astype(int), Lk.astype(int)


def calc_BNkM(AMNk):
    # TODO implement simplier method from book
    return np.linalg.inv(AMNk)


def starting_vector(cf:NpCanonicalForm) -> np.array:
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
            return xN
    print("error, initial vector not found!")
    return np.zeros(cf.n)


def simplex_step(cf: NpCanonicalForm, xkN: np.array) -> (np.array, bool):
    """
    Совершает один шаг алгоритма симплекс-метода. Обозначения и процедура взяты из пособия
    Петухов и др., стр. 88
    :param cf: параметры задачи, поставленной в канонической форме
    :param xkN: начальное приближение - некий опорный вектор к множеству, заданному `cf`
    :return:
        1. следующее приближение - тоже опорный вектор, причём с ним значение целевой функции не возрастает
        2. булеву переменную, равную `true`, если итерирование нужно прекратить и `false` иначе
    """
    AMN, cN = cf.A, cf.c
    _, _, Nk0, Nk_plus = split_xkN(xkN) # Nk0 - индексы нулевых компонент xk, Nk+ - положительных
    binomGrid = binomial_grid(len(Nk0), cf.m - len(Nk_plus))  # вспомог. структура для построения комбинаций столбцов, присоединяемых к A[M,Nk+]
    for binom_idx in range(binomGrid[-1, -1]):  # итерируемся по комбинациям векторов, присоединяемых к A[M,Nk+]
        AMNk, Nk, Lk = new_AMNk(AMN, xkN, binomGrid, binom_idx)  # дополняем Nk+ до Nk так, что A[M, Nk] квадратная
        if np.linalg.det(AMNk) == 0:  # если определитель построенной квадратной матрицы 0, пропускаем комбинацию
            continue
        BNkM = calc_BNkM(AMNk)  # вычисление матрицы, обратной к A[M, Nk]
        cNk = np.array([cN[i] for i in Nk])
        ykM = np.matmul(BNkM.T, cNk)
        dkN = cN - np.matmul(AMN.T, ykM)
        dkLk = np.array([dkN[int(i)] for i in Lk])
        if np.min(dkLk) >= 0:  # xkN уже является оптимальным вектором
            print("solution found!")
            return xkN, True
        jk = list(filter(lambda j: dkLk[j] < 0, range(len(Lk))))[0]  # индекс первой негативной компоненты в dkLk
        xkNk0, xkNk_plus, Nk0, Nk_plus = split_xkN(xkN)
        ukNk = np.matmul(BNkM, AMN[:,jk])
        if np.max(ukNk) <= 0:  # целевая функция не ограничена снизу
            print("solution does not exist")
            return np.array([np.inf for _ in range(cf.n)]), True
        if len(Nk_plus) == len(Nk) or max([ukNk[i] for i in filter(lambda j: Nk[j] not in Nk_plus, range(len(Nk)))]) < 0:
            ukN = [ukNk[list(Nk).index(i)] if i in Nk else 0 for i in range(cf.n)]
            ukN[jk] = -1
            theta_k = min([xkN[i]/ukN[i] for i in filter(lambda j: ukN[j] > 0, Nk)])
            return xkN - np.multiply(theta_k, ukN), False
    return xkN, False
