"""
command line launcher for gradient descent methods
usage: from containing folder:
python -m src.grad_cli --out src/out.txt --x0 0. --y0 0. --alg newton
"""
import getopt
import sys
from typing import List
import numpy as np

from src import GradientDescent, Newton


def run_steepest_descent(x0: float, y0: float) -> List[np.array]:
    return GradientDescent.run_gradient_descent(x0=np.array([x0, y0]))


def run_newton(x0: float, y0: float) -> List[np.array]:
    return Newton.newton(np.array([x0, y0]))


def parse(argv):
    try:
        opts, args = getopt.getopt(argv, "", ["out=", "x0=", "y0=", "alg="])
    except getopt.GetoptError as err:
        print(err)
        return
    kwargs = {}
    for opt in opts:
        kwargs[opt[0]] = opt[1]
    x0 = float(kwargs["--x0"]) if "--x0" in kwargs else 0.
    y0 = float(kwargs["--y0"]) if "--y0" in kwargs else 0.
    output_filename = kwargs["--out"] if "--out" in kwargs else "steps.txt"
    if "--alg" in kwargs and kwargs["--alg"] == "newton":
        steps = run_newton(x0, y0)
    else:
        steps = run_steepest_descent(x0, y0)
    with open(output_filename, "w+") as file_out:
        for x in steps:
            line = str(x[0]) + " " + str(x[1]) + "\n"
            file_out.write(line)


if __name__ == "__main__":
    parse(sys.argv[1:])
