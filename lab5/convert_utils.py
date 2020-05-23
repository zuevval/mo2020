import numpy as np


class LinearProblem:
    """
    c^T*x->min
    x in S:={x |A*x<=b}
    (a particular case of general linear programming problem)
    """
    def __init__(self, A: np.ndarray, b: np.array, c:np.array):
        if len(A) != len(b):
            raise Exception("dimensions mismatch")
        for line in A:
            if len(line) != len(c):
                raise Exception("dimensions mismatch")
        self.c, self.A, self.b = c, A, b


class CanonicalForm:
    """
    c^T*x->min
    x in S:={x|A*x=b, x >= 0}
    """
    def __init__(self, A:np.ndarray, b:np.array, c:np.array):
        self.n = len(c)
        self.m = len(b)
        dimensions_matched = len(A) == self.m
        for row in A:
            dimensions_matched &= len(row) == self.n
        if not dimensions_matched:
            raise Exception("dimensions mismatch")
        self.A, self.b, self.c = A, b, c


def dual_problem(lp: LinearProblem) -> CanonicalForm:
    A_dual = - lp.A.T
    b_dual = lp.c
    c_dual = lp.b
    return CanonicalForm(A=A_dual, b=b_dual, c=c_dual)


def back_to_primal(Nk: np.array, primal: LinearProblem) -> np.array:
    return primal.b[Nk].dot(np.linalg.inv(primal.A[Nk]))


def example():
    lp = LinearProblem(
        A=np.array([[1, 2, 3], [4, 5, 6], [3, 2, 1], [6, 4, 5]]),
        b=np.array([0, 1, 2, 3]),
        c=np.array([5, 3, 1])
    )
    print(lp.A)
    dcf = dual_problem(lp)
    print(dcf.A)


if __name__ == '__main__':
    example()
