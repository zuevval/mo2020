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
    x in S:={x|A*x=b}
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


def lp_to_canonical(lp: LinearProblem, func, MatrixSigns, VariablesSigns) -> CanonicalForm:
    cf = cn.Convert(lp.A,lp.b,lp.c, func, MatrixSigns, VariablesSigns)
    return cf

def dual_problem(lp: LinearProblem, func, MatrixSigns, VariablesSigns) -> CanonicalForm:
    dual_cf = dl.CreateDual(lp.A,lp.b,lp.c, func, MatrixSigns, VariablesSigns)
    return dual_cf

def back_to_primal(x: np.array, primal: CanonicalForm, dual: CanonicalForm) -> np.array:
    # TODO given a solution of dual problem, restore solution of primal
    pass

if __name__ == '__main__':
    A,b,c,func,MatrixSigns,VariablesSigns = et.ReadFile('problem.txt')
    A = np.array(A)
    b = np.array(b)
    c = np.array(c)
    lp = LinearProblem(
        A,
        b,
        c
    )
    cf = lp_to_canonical(lp, func, MatrixSigns, VariablesSigns)
    dcf = dual_problem(lp, func, MatrixSigns, VariablesSigns)