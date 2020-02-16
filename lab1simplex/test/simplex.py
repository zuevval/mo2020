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


def prepare_common_data():
    A, b, svecs = prepare_data_from_lectures()
    c = [1, 2, 3, 4]
    ncf = utils.NpCanonicalForm(utils.CanonicalForm(A, b, c))
    xNk = np.array(svecs[0])
    n_tests = 10
    return ncf, svecs, xNk, n_tests


def test_cost_function_doesnt_increase():
    """
    ensures that cost function (c^T*x) doesn't increase over iterations
    """
    ncf, _, xNk, n_tests = prepare_common_data()
    for _ in range(n_tests):
        xNk_minus1, xNk = xNk, simplex.simplex_step(ncf, xNk)
        assert np.dot(xNk_minus1, ncf.c) <= np.dot(xNk, ncf.c)


def test_if_vector_in_s():
    """
    checks whether a result of a simplex_step(...) is in S={x|Ax=b}
    """
    ncf, _, xNk, n_tests = prepare_common_data()
    abs_error = 1e-3
    for _ in range(n_tests):
        xNk = simplex.simplex_step(ncf, xNk)
        assert np.sum(np.abs(np.matmul(ncf.A, xNk) - ncf.b)) < abs_error

def test_if_vector_support():
    """
    makes sure that a result of a simplex_step(...) belongs to the set of support vectors
    """
    ncf, svecs, xNk, n_tests = prepare_common_data()
    abs_error = 1e-3
    for _ in range(n_tests):
        xNk = simplex.simplex_step(ncf, xNk)
        dists = np.array([np.sum(np.square(np.array(sv) - xNk)) for sv in svecs])
        assert np.min(dists) < abs_error  # assert that exists a support vector close to xNk


def run_all_simplex_tests():
    test_cost_function_doesnt_increase()
    test_if_vector_in_s()
    test_if_vector_support()


if __name__ == "__main__":
    run_all_simplex_tests()