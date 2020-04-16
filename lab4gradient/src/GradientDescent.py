from src.GoldenSelection import golden_ratio
from numpy import linalg, dot, array
from typing import Callable, NewType, List
from math import cos, sin
import sys

Point = NewType('Point', array)
Vector = NewType('Vector', array)


def sample_function(
        x: Point
) -> float:
    """
        Vector function used for testing.
        in: x: np.array -  point in R^n
        out: test_function(x): float - function at given x
    """
    return pow(x[0], 2) + pow(x[1], 2) + cos(x[0] + 3 * x[1]) - x[0] + 2 * x[1]


def sample_function_gradient(
        x: Point
) -> Vector:
    """
        Gradient of test_function above.
        in: x: np.array - point in R^n
        out: test_function_gradient(x): np.array - vector of gradient at given x
    """
    return Vector([-sin(x[0] + 3 * x[1]) + 2 * x[0] - 1, -3 * sin(x[0] + 3 * x[1]) + 2 * x[1] + 2])


def gradient_descent(
        function: Callable[[Point], float],
        gradient: Callable[[Point], Vector],
        x_0: Point,
        eps: float
) -> Point:
    """
        Minimizes n-dimensional function.
        in: function: Callable[[np.array], float] - function to minimize
            gradient: Callable[[np.array], np.array] - gradient of function f
            x0: np.array - starting point
            eps: float - stopping criterion
        out: x_min: np.array - argument such that ||nabla(f(x_min))||^2 < eps
    """
    x_k = x_0
    while True:
        yield x_k
        grad_x_k = gradient(x_k)
        if pow(linalg.norm(grad_x_k), 2) < eps:
            return x_k
        alpha_k = golden_ratio(lambda alpha: function(x_k - dot(alpha, grad_x_k)), (0, 10), eps)
        x_k = x_k - dot(alpha_k, grad_x_k)


def run_gradient_descent(f: Callable[[array], float] = sample_function,
                         grad: Callable[[array], array] = sample_function_gradient,
                         x0: array = array([0., 0.]), eps: float = 0.1,
                         output_filename: str = "steps.txt") -> List[array]:
    """
    runs steepest descent minimization algorithm and logs position on every
    iteration to output file
    :param f: function to be minimized, f: R^n->R
    :param grad: gradient of `f`, grad: R^n->R^n
    :param x0: starting point, x0 in R^n
    :param eps: exit criteria (when norm(grad) < eps - stop)
    :param output_filename: path to desired output file including extension
    :return: array of points visited while running algorithm, each in R^n
    """
    steps = list(gradient_descent(f, grad, x0, eps))
    with open(output_filename, "w+") as file_out:
        for x in steps:
            line = str(x[0]) + " " + str(x[1]) + "\n"
            file_out.write(line)
    return steps


if __name__ == "__main__":
    if len(sys.argv) > 1:  # if output filename passed as command line argument
        out_filename = str(sys.argv[1])

        def f_simple(x: array) -> float:
            return pow(x[0], 2) + pow(x[1], 2) - x[0] + 2 * x[1]

        def f_simple_grad(x: array) -> array:
            return array([2 * x[0] - 1, 2 * x[1] + 2])

        run_gradient_descent(f_simple, f_simple_grad,
                             array([-0.5, -0.5]), output_filename=out_filename)
    else:
        res = run_gradient_descent(sample_function, sample_function_gradient,
                                   array([0, 0]),
                                   output_filename="real_function_steps.txt")
        for point in res:
            print(point)
