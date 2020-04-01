from math import sin, cos, pi
def function(x, counter):
    counter = counter + 1
    return sin(pi*x) + cos(pi*x), counter


def FibonacciSeq(n):
    if n in (1, 2):
        return 1
    return FibonacciSeq(n - 1) + FibonacciSeq(n - 2)


def FibonacciMethod(a_k, b_k, n, k, lambda_k, mu_k, delta, prev_function, counter, eps):
    if k <= n - 2:
        if lambda_k==float('inf'):
            lambda_k = a_k + FibonacciSeq(n - k - 1) / FibonacciSeq(n - k + 1) * (b_k - a_k)
            if k != 0:
                f1, counter = function(lambda_k, counter)
                f2 = prev_function
        if mu_k==float('inf'):
            mu_k = a_k + FibonacciSeq(n - k) / FibonacciSeq(n - k + 1) * (b_k - a_k)
            if k != 0:
                f1 = prev_function
                f2, counter = function(mu_k, counter)
        if k == 0:
            f1, counter = function(lambda_k, counter)
            f2, counter = function(mu_k, counter)
    else:
        print("l =", b_k - a_k)
        print("x_min = ", (a_k + b_k) / 2, "\nFunction addressed", counter + 1, "times")
        return

    
    if f1 > f2:
        prev_function = f2
        if k < n - 2:
            FibonacciMethod(lambda_k, b_k, n, k + 1, mu_k, float('inf'), delta, prev_function, counter, eps)
        else:
            FibonacciMethod(lambda_k, b_k, n, k + 1, mu_k, lambda_k, delta, prev_function, counter, eps)
    else:
        prev_function = f1
        if k < n - 2:
            FibonacciMethod(a_k, mu_k, n, k + 1, float('inf'), lambda_k, delta, prev_function, counter, eps)
        else:
            FibonacciMethod(lambda_k, b_k, n, k + 1, mu_k, lambda_k, delta, prev_function, counter, eps)


if __name__=="__main__":
    counter = 0  # count how many times function was addressed
    a = 0
    b = 2
    delta = 0.01
    eps = 0.001
    N = (b - a) / eps
    n = 1
    while FibonacciSeq(n) < N:
        n = n + 1

    print("Fibonacci method:")
    FibonacciMethod(a, b, n, 0, float('inf'), float('inf'), delta, float('inf'), counter, eps)
