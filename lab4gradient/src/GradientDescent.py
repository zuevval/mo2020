from src.GoldenSelection import golden_ratio
from numpy import linalg, dot, array
from typing import Callable, NewType
from math import cos, sin

Point = NewType('Point', array)
Vector = NewType('Vector', array)


def test_function(
        x: Point
) -> float:
    """
        Vector function used for testing.
        in: x: np.array -  point in R^n
        out: test_function(x): float - function at given x
    """
    return pow(x[0], 2) + pow(x[1], 2) + cos(x[0] + 3 * x[1]) - x[0] + 2 * x[1]


def test_function_gradient(
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
        grad_x_k = gradient(x_k)
        if pow(linalg.norm(grad_x_k), 2) < eps:
            return x_k
        alpha_k = golden_ratio(lambda alpha: function(x_k - dot(alpha, grad_x_k)), (0, 10), eps)
        x_k = x_k - dot(alpha_k, grad_x_k)


if __name__ == "__main__":
    x_0 = array([0, 0])
    eps = 0.1
    print(gradient_descent(test_function, test_function_gradient, x_0, eps))