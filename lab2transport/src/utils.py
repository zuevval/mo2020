import numpy as np


class TransportTable:
    def __init__(self, u: np.array, v: np.array, c: np.ndarray, x0: np.ndarray, alpha: np.ndarray):
        self.m, self.n = len(u), len(v)
        assert c.shape == x0.shape == alpha.shape == (self.m, self.n)
        assert len([xij for xij in x0.flatten() if xij is not None]) == self.m + self.n - 1
        self.c, self.x, self.alpha, self.u, self.v = c, x0, alpha, u, v
