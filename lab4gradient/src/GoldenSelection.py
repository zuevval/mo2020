import math
from typing import Callable, NewType, Tuple

Interval = NewType('Interval', Tuple[int, int])


def counting(func: Callable):
    """
        Decorator that keeps track of the number of times a function is called.
        in: func: Callable - function to wrap
        out: ncalls: int - number of calls as a attribute of f
    """

    def wrapper(*args, **kwargs):
        wrapper.ncalls += 1
        # print("{0} called: {1}x".format(func.__name__, wrapper.ncalls))
        return func(*args, **kwargs)

    wrapper.ncalls = 0
    return wrapper


def golden_ratio(f: Callable[[float], float], ab: Interval, eps: float) -> Interval:
    """
        Golden section search.

        Given a function f with a single local minimum in
        the interval [a,b], goldenRatio returns a subset
        interval of length <= eps that contains the minimum.

        in: f: Callable[[float], float] - function to minimize
            ab: Tuple[float, float] - interval with a single local minimum of f
            eps: float - desired length of uncertainty interval
        out: cd: Tuple[float, float] - uncertainty interval
    """
    invphi2 = 0.5 * (3 - math.sqrt(5))  # 1/phi^2

    (a_k, b_k) = (min(ab), max(ab))
    interval_length = b_k - a_k
    if interval_length <= eps:
        return a_k, b_k

    lambda_k = a_k + invphi2 * interval_length
    mu_k = b_k - invphi2 * interval_length
    f_lambda_k = f(lambda_k)
    f_mu_k = f(mu_k)
    assert (mu_k - a_k == b_k - lambda_k), "Length of new uncertainty interval depends on the comparison interval!"

    while interval_length > eps:
        if f_lambda_k < f_mu_k:
            b_k = mu_k
            mu_k = lambda_k
            f_mu_k = f_lambda_k
            interval_length = b_k - a_k
            lambda_k = a_k + invphi2 * interval_length
            f_lambda_k = f(lambda_k)
        else:
            a_k = lambda_k
            lambda_k = mu_k
            f_lambda_k = f_mu_k
            interval_length = b_k - a_k
            mu_k = b_k - invphi2 * interval_length
            f_mu_k = f(mu_k)

    #uncertainty interval
    start, end = (a_k, mu_k) if f_lambda_k < f_mu_k else (lambda_k, b_k)
    return (end - start)*0.5
    #return (a_k, mu_k) if f_lambda_k < f_mu_k else (lambda_k, b_k)


