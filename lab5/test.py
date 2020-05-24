from convert_utils import dual_problem, back_to_primal, LinearProblem
from cutting_plane import cutting_plane_alg
import numpy as np


def prepare_simple_data() -> LinearProblem:
    return LinearProblem(A=np.array([[-2],
                                   [-1],
                                   [-2]]),
                         b=np.array([-1, -1, -1]),
                         c=np.array([2]))


def test_convert_utils():
    lp = prepare_simple_data()
    dual = dual_problem(lp)
    assert (dual.A == np.array([[2, 1, 2]])).all()
    assert (dual.b == lp.c).all()
    assert (dual.c == lp.b).all()
    dual_solution = np.array([0, 2, 0])
    Nk_dual = np.array([1])
    x_primal = back_to_primal(Nk_dual, lp)
    assert (x_primal == np.array([1])).all()


def test_cutting_plane():
    lp = prepare_simple_data()
    x = cutting_plane_alg(lp=lp,
                          phi=lambda x: np.dot(lp.c, x),
                          phi_subgrad=lambda _: lp.c)
    assert (x == np.array([1])).all()

