from typing import List

class CanonicalForm:
    def __init__(self, A:List[List[float]], b:List[float], c:List[float]):
        self.n = len(c)
        self.m = len(b)
        dimensions_matched = len(A) == self.m
        for row in A:
            dimensions_matched &= len(row) == self.n
        if not dimensions_matched:
            raise Exception("dimensions mismatch")
        self.A, self.b, self.c = A, b, c
