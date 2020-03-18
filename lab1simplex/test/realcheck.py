"""
Checks performed with real parameters
"""

from src import enter, canon, utils, simplex, bruteforce
import numpy as np

def run_primal_problem():
    cf = utils.NpCanonicalForm(canon.Convert(*enter.ReadFile("../src/problem.txt")))
    x0 = simplex.starting_vector(cf)
    assert np.max(np.abs(np.matmul(cf.A, x0) - cf.b)) < 1e-3
    print(np.matmul(x0, cf.c))
    res = simplex.simplex_method(cf, x0)
    # assert np.max(np.abs(np.matmul(cf.A, res.x) - cf.b)) < 1e-3
    f = np.matmul(res.x, cf.c)
    print("cost function at vector found by simplex: " + str(f))

    x_brute = bruteforce(cf)
    assert np.max(np.abs(np.matmul(cf.A, x_brute) - cf.b)) < 1e-3
    f_b = np.matmul(x_brute, cf.c)
    print("cost function at vector found by brute force: " + str(f_b) + "\n")


def check_with_automatically_converted_canonical():
    c = [8, 8, -8, -4, 4, -2, 2, 0, 0]
    A = [[1, -2, 2, 2, -2, 0, 0, 1, 0],
         [1, 2, -2, 1, -1, 1, -1, 0, 0],
         [2, 1, -1, -4, 4, 1, -1, 0, 0],
         [-1, 4, -4, 0, 0, -2, 2, 0, 1]]
    b = [6, 24, 30, -6]
    cf = utils.NpCanonicalForm(utils.CanonicalForm(A, b, c))

    x_online = np.array([0, 0, 7 / 2, 0, 1 / 2, 63 / 2, 0, 0, 71])
    print(np.matmul(x_online, cf.c))
    print(str(np.matmul(cf.A, x_online) - cf.b))

    print("\n--- running test with pre-defined canonical form ---")
    print("---- running simplex algorithm ----")
    res = simplex.simplex_method(cf, simplex.starting_vector(cf))
    print(res.x)
    f = np.matmul(res.x, cf.c)
    print("cost function at vector found by simplex: " + str(f))

    print("\n---- running brute force ----")
    x_b = bruteforce(cf)
    print(x_b)
    assert np.max(np.abs(np.matmul(cf.A, x_b) - cf.b)) < 1e-3
    f_b = np.matmul(x_b, cf.c)
    print("cost function at vector found by brute force: " + str(f_b))


if __name__ == "__main__":
    run_primal_problem()
    check_with_automatically_converted_canonical()

