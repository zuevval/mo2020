"""
Checks performed with real parameters
"""

from src import enter, canon, utils, simplex, bruteforce
import numpy as np

def run_primal_problem():
    cf = utils.NpCanonicalForm(canon.Convert(*enter.ReadFile("../src/problem.txt")))
    x0 = simplex.starting_vector(cf)
    res = simplex.simplex_method(cf, x0)
    assert np.max(np.abs(np.matmul(cf.A, res.x) - cf.b)) < 1e-3
    f = np.matmul(res.x, cf.c)
    print(f)

    x_brute = bruteforce(cf)
    assert np.max(np.abs(np.matmul(cf.A, res.x) - cf.b)) < 1e-3
    f_b = np.matmul(x_brute, cf.c)
    print(f_b)



if __name__ == "__main__":
    run_primal_problem()

