import math


def counting(func):
    """ Decorator that keeps track of the number of times a function is called """

    def wrapper(*args, **kwargs):
        wrapper.ncalls += 1
        # print("{0} called: {1}x".format(func.__name__, wrapper.ncalls))
        return func(*args, **kwargs)

    wrapper.ncalls = 0
    return wrapper


def goldenRatio(f, a, b, eps):
    invphi2 = 0.5 * (3 - math.sqrt(5))  # 1/phi^2

    (a_k, b_k) = (min(a, b), max(a, b))
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

    return (a_k, mu_k) if f_lambda_k < f_mu_k else (lambda_k, b_k)


@counting
def testFunction(x):
    return pow(x, 2) + 2 * x


if __name__ == "__main__":
    a = -3
    b = 5
    eps = 0.00001
    print("Uncertainty interval is " + str(goldenRatio(testFunction, a, b, eps)))
    print("Function addressed " + str(testFunction.ncalls) + " times.")
