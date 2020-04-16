from typing import List

from math import cos, sin, sqrt
import numpy as np

def f(x1,x2):
    return pow(x1,2)+pow(x2,2)+cos(x1+3*x2)-x1+2*x2

def F(x1,x2):
    return 2*x1-sin(x1+3*x2)-1, 2*x2-3*sin(x1+3*x2)+2

def HesseCoef(x1,x2):
    F11 = 2 - cos(x1+3*x2)
    F12 = -3*cos(x1+3*x2)
    F21 = -3*cos(x1+3*x2)
    F22 = 2 - 9*cos(x1+3*x2)
    return F11,F12,F21,F22

def det(F11,F12,F21,F22):
    return F11*F22-F21*F12


def newton(x_start:np.array = np.array([1., 1.])) -> List[np.array]:
    x1 = x_start[0]
    x2 = x_start[1]
    res = [x_start]
    eps = 0.1
    F1, F2 = F(x1,x2)
    while sqrt(pow(F2,2)+pow(F1,2)) > eps:
        F11,F12,F21,F22 = HesseCoef(x1, x2)
        detH = det(F11,F12,F21,F22)
        F1, F2 = F(x1,x2)
        x1 = x1 - 1/detH*(F22*F1-F21*F2)
        x2 = x2 - 1/detH*(F11*F2-F12*F1)
        res.append(np.array([x1, x2]))

    print("x1 = ", x1,"x2 = ", x2, "f(x1,x2) =", f(x1,x2))
    return res

if __name__ == "__main__":
    newton()