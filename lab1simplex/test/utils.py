from src.utils import *


def test_binomial_grid():
    n, k = 5, 2
    bg = binomial_grid(n, k)
    assert len(bg) == n - k + 1
    assert len(bg[0]) == k + 1
    assert bg[n - k, k] == 10
    for row in range(n - k + 1):
        assert bg[row, 0] == 1
    for col in range(k + 1):
        assert bg[0, col] == 1


def test_subset_by_index():
    n, k = 5, 3
    set = np.array(list(range(n)))
    bg = binomial_grid(n, k)
    subsets = []
    for idx in range(bg[-1, -1]):
        subsets.append(subset_by_index(set, bg, idx))
    print(subsets)


def run_all_utils_tests():
    test_binomial_grid()
    test_subset_by_index()


if __name__ == "__main__":
    run_all_utils_tests()