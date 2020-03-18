import numpy as np
from typing import Tuple

def closeForm(matrix:np.ndarray, stock:np.array, need:np.array)-> Tuple[np.ndarray, np.array, np.array]:
    """Brings the transportation problem to closed-form by adding fake provider or consumer"""
    overallNeed = np.sum(need)
    overallStock = np.sum(stock)
    rows, cols = matrix.shape
    delta = overallStock - overallNeed
    if delta < 0:
        fakeRow = np.zeros(shape=(1, cols))
        matrix = np.append(matrix, fakeRow, axis=0)
        stock = np.append(stock, -delta)
    if delta > 0:
        fakeColumn = np.zeros(shape=(rows, 1))
        matrix = np.append(matrix, fakeColumn, axis=1)
        need = np.append(need, delta)
    return matrix, stock, need