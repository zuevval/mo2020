from src import utils, simplex
import numpy as np

# ----------- Unit Tests --------------

def prepare_data_from_lectures():
    A = [[5, 5, 6, 4],
         [4, 5, -5, -4]]
    b = [5, 4]
    support_vectors = [[1, 0, 0, 0],
                       [0, 49/55, 1/11, 0],
                       [0, 9/10, 0, 1/8]]
    return A, b, support_vectors


def prepare_common_data(svecs_idx:int):
    A, b, svecs = prepare_data_from_lectures()
    c = [1, 2, 3, 4]
    ncf = utils.NpCanonicalForm(utils.CanonicalForm(A, b, c))
    xkN = np.array(svecs[svecs_idx])
    n_tests = 10
    return ncf, svecs, xkN, n_tests


def test_split_xkN():
    for xkN in [[0],
                [0, 1, 2, 3],
                [1, 2, 3, 4],
                [1, 0, 1, 0],
                [0, 0, 3.5, 1.1, 0, 1, 2, 0, 0, 0.01]]:
        xkNk0, xkNk_plus, Nk0, Nk_plus = simplex.split_xkN(xkN)
        for i in range(len(xkN)):
            if xkN[i] > 0:
                assert xkN[i] in xkNk_plus
                assert xkN[i] not in xkNk0
                assert i in Nk_plus
                assert i not in Nk0
            else:
                assert xkN[i] in xkNk0
                assert xkN[i] not in xkNk_plus
                assert i in Nk0
                assert i not in Nk_plus


def test_new_AMNk_simple():
    abs_err = 1e-3
    ncf, _, _, _ = prepare_common_data(0)
    AMN = ncf.A
    x = [0, 1, 0, 1]
    AMNk, _, _ = simplex.new_AMNk(AMN, x, utils.binomial_grid(2, ncf.m - 2), 0)
    for i, j in zip([0, 1], [1, 3]):
        for row in range(len(AMN)):
            assert AMNk[row, i] == AMN[row, j]
    assert abs(np.linalg.det(AMNk)) > abs_err

    x = [1, 0, 0, 0]
    AMNk, _, _ = simplex.new_AMNk(AMN, x, utils.binomial_grid(3, ncf.m - 1), 0)
    for i, j in zip([0, 1], [0, 3]):
        for row in range(len(AMN)):
            assert AMNk[row, i] == AMN[row, j]
    assert abs(np.linalg.det(AMNk)) > abs_err
    AMNk, _, _ = simplex.new_AMNk(AMN, x, utils.binomial_grid(3, ncf.m - 1), 2)
    for i, j in zip([0, 1], [0, 1]):
        for row in range(len(AMN)):
            assert AMNk[row, i] == AMN[row, j]
    assert abs(np.linalg.det(AMNk)) > abs_err

    x = [0, 0, 1, 0]
    AMNk, _, _ = simplex.new_AMNk(AMN, x, utils.binomial_grid(3, ncf.m - 1), 2)
    for i, j in zip([0, 1], [0, 2]):
        for row in range(len(AMN)):
            assert AMNk[row, i] == AMN[row, j]
    assert abs(np.linalg.det(AMNk)) > abs_err


def test_calc_BNkM():
    abs_err = 1e-3
    ncf, _, _, _ = prepare_common_data(0)
    A1 = ncf.A[:, [0, 1]]
    A2 = ncf.A[:, [1, 2]]
    A3 = ncf.A[:, [2, 3]]
    E = np.eye(2)
    for Ai in [A1, A2, A3]:
        assert abs(np.max(np.matmul(Ai, simplex.calc_BNkM(Ai)) - E)) < abs_err


def test_starting_vector():
    cf = prepare_common_data(0)[0]
    print("test not implemented")

# --------- Integration Tests ------------

def test_cost_function_doesnt_increase():
    """
    ensures that cost function (c^T*x) doesn't increase over iterations
    """
    for svecs_idx in range(3):
        ncf, _, xkN, n_tests = prepare_common_data(svecs_idx)
        for _ in range(n_tests):
            xk_minus1_N, (xkN, _, stopIteration) = xkN, simplex.simplex_step(ncf, xkN)
            if stopIteration:
                break
            assert np.dot(xk_minus1_N, ncf.c) >= np.dot(xkN, ncf.c)


def test_if_vector_in_s():
    """
    checks whether a result of a simplex_step(...) is in S={x|Ax=b, x >= 0}
    """
    for svecs_idx in range(3):
        ncf, _, xkN, n_tests = prepare_common_data(svecs_idx)
        abs_error = 1e-3
        for _ in range(n_tests):
            xkN, _, stopIteration = simplex.simplex_step(ncf, xkN)
            if stopIteration:
                break
            assert np.sum(np.abs(np.matmul(ncf.A, xkN) - ncf.b)) < abs_error
            assert min(xkN) >= 0  # assert there are no negative components in xkN


def test_if_vector_support():
    """
    makes sure that a result of a simplex_step(...) belongs to the set of support vectors
    """
    for svecs_idx in range(3):
        ncf, svecs, xkN, n_tests = prepare_common_data(svecs_idx)
        abs_error = 1e-3
        for _ in range(n_tests):
            xkN, _, stopIteration = simplex.simplex_step(ncf, xkN)
            if stopIteration:
                break
            dists = np.array([np.sum(np.square(np.array(sv) - xkN)) for sv in svecs])
            assert np.min(dists) < abs_error  # assert that exists a support vector close to xkN


def test_simplex_method():
    abs_error = 1e-3
    for svecs_idx in range(3):
        ncf, svecs, xkN, _ = prepare_common_data(svecs_idx)
        res = simplex.simplex_method(ncf, xkN)
        assert res.is_solution
        dists = np.array([np.sum(np.square(np.array(sv) - res.x)) for sv in svecs])
        assert np.min(dists) < abs_error  # assert that exists a support vector close to xkN

def run_all_simplex_tests():
    # unit tests
    test_split_xkN()
    test_new_AMNk_simple()
    test_calc_BNkM()
    test_starting_vector()
    # integration tests
    test_cost_function_doesnt_increase()
    test_if_vector_in_s()
    test_if_vector_support()
    test_simplex_method()


if __name__ == "__main__":
    run_all_simplex_tests()
