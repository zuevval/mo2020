import numpy as np
import canon as cn
import enter as et
import dual as dl


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
    A1 = np.transpose(lp.A)
    c1 = lp.b
    b1 = -lp.c
    return CanonicalForm(A=A1, b=b1, c=c1)


def back_to_primal(x: np.array, inv_AMNk: np.array, Nk:np.array,
                   primal: LinearProblem, dual: CanonicalForm) -> np.array:
    # TODO remove unused variables
    """
    `dual` is a canonical form of linear problem
    (c[N],x[N])->min, x in S:={x|A[M,N]x[N]<=b[M]}
    and `primal` is a form with
    (c1[M]=-b[M], y[M])->min, y in S1:={y|A1[N,M] = b1[N]}
    A1 = A^T - transposed A, b1 = c

    `x` is a solution of `dual`, `Nk` - its basis, inv_AMNk = A[M,Nk]^-1
    """
    # variable names below: see report 1, "primal problem solution restoration"
    c = dual.b
    return np.matmul(c, inv_AMNk)


if __name__ == '__main__':
    A, b, c, _, _, _ = et.ReadFile('problem.txt')
    A = np.array(A)
    b = np.array(b)
    c = np.array(c)
    lp = LinearProblem(
        A,
        b,
        c
    )
    print(A)
    dcf = dual_problem(lp)
    print(dcf.A)
