from src import utils, simplex
import numpy as np


def prepare_data_from_lectures():
    A = [[5, 5, 6, 4],
         [4, 5, -5, -4]]
    b = [5, 4]
    support_vectors = [[1, 0, 0, 0],
                       [0, 49/55, 1/11, 0],
                       [0, 9/10, 0, 1/8]]
    return A, b, support_vectors


def test_cost_function_not_decreases():
    A, b, svecs = prepare_data_from_lectures()
    c = [1, 2, 3, 4]
    ncf = utils.NpCanonicalForm(utils.CanonicalForm(A, b, c))
    xNk = svecs[0]
    n_tests = 10
    for _ in range(n_tests):
        xNk_minus1, xNk = xNk, simplex.simplex_step(ncf, xNk)
        assert np.dot(xNk_minus1, ncf.c) <= np.dot(xNk, ncf.c)


def run_all_simplex_tests():
    test_cost_function_not_decreases()


if __name__ == "__main__":
    run_all_simplex_tests()