from typing import List
import numpy as np


class CanonicalForm:
    """
    Класс, хранящий входные данные канонической формы задачи линейного программирования.
    Обозначения переменных: пособие Петухова и др., стр. 80, ф.-ла (4.25)
    """
    def __init__(self, A:List[List[float]], b:List[float], c:List[float]):
        self.n = len(c)
        self.m = len(b)
        dimensions_matched = len(A) == self.m
        for row in A:
            dimensions_matched &= len(row) == self.n
        if not dimensions_matched:
            raise Exception("dimensions mismatch")
        self.A, self.b, self.c = A, b, c


class NpCanonicalForm:
    """
    Представление канонической формы в виде массивов NumPy (в отличие от списков в `CanonicalForm`)
    """
    def __init__(self, cf:CanonicalForm):
        self.n, self.m = cf.n, cf.m
        self.A = np.array([np.array(row) for row in cf.A])
        self.b = np.array(cf.b)
        self.c = np.array(cf.c)


def binomial_grid(n: int, k: int) -> np.ndarray:
    """
    Конструирует прямоугольную матрицу `n-k+1`*`k+1` - повёрнутый набок кусок треугольника Паскаля
    с углами в (C[0,0], C[`n`,`k`])
    :param n: основание биномиального коэффициента C[`n`,`k`]
    :param k: номер биномиального коэффициента C[`n`,`k`]
    :return: матрицу с биномиальными коэффициентами
    """
    res = [[1 for _ in range(k+1)]]
    for row in range(n-k):
        res.append([1])
        for col in range(1,k+1):
            res[-1].append(res[-1][-1] + res[-2][col])  # С[row,col]=C[row,col-1]+C[row-1,col]
    return np.array([np.array(row) for row in res])


def subset_by_index(set:np.array, binom_table:np.ndarray, binom_index:int) -> np.array:
    """
    Возвращает подмножество заданного массива длины `n` по номеру от 0 до C[n,k],
    где C[n,k] - число, стоящее в правом нижнем углу `binom_table`, k < n.
    Каждому номеру `binom_index` алгоритм ставит в соответствие подмножество, для его вычисления используется таблица.
    :param set: исходный массив элементов
    :param binom_table: биномиальная таблица, созданная функцией `binomial_grid`
    :param binom_index: номер набора
    :return: подмножество (элементы в том же порядке, что и в `set`)
    """
    res, row, col = [], len(binom_table) - 1, len(binom_table[-1]) - 1
    while col > 0:
        if binom_index < binom_table[row, col - 1]:
            res.append(set[row + col - 1])
            col -= 1
        else:
            binom_index -= binom_table[row, col - 1]
            row -= 1
    res = np.array(res)
    res.sort()
    return res.astype(int)


if __name__ == "__main__":
    inputData = CanonicalForm([[1, 2, 3], [4, 5, 6]], [7, 8], [9, 10, 11])
    npData = NpCanonicalForm(inputData)
    print(npData.A)
    print(npData.b)
    print(npData.c)
    print(npData.n)

    print(binomial_grid(5, 2))
