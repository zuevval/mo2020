from cutting_plane import cutting_plane_alg as cut_alg
from convert_utils import LinearProblem
import numpy as np
import logging
import enter
import sys

if __name__ == "__main__":
    logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
    A, b, c, _, _, _ = enter.ReadFile('problem.txt')
    lp = LinearProblem(np.array(A), np.array(b), np.array(c))
    # TODO launch cutting plane algorithm
