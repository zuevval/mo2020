from typing import List
import numpy as np


class CanonicalForm:
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
    def __init__(self, cf:CanonicalForm):
        self.n, self.m = cf.n, cf.m
        self.A = np.array([np.array(row) for row in cf.A])
        self.b = np.array(cf.b)
        self.c = np.array(cf.c)


if __name__ == "__main__":
    inputData = CanonicalForm([[1, 2, 3], [4, 5, 6]], [7, 8], [9, 10, 11])
    npData = NpCanonicalForm(inputData)
    print(npData.A)
    print(npData.b)
    print(npData.c)
    print(npData.n)