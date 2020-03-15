import numpy as np
import sys
sys.path.append('.')  # now environment can see package `src`
from src import NW
from src import Min

if __name__ == "__main__":
    matrix1 = np.array([[5,3,1],[3,2,4],[4,1,2]], int)
    stock1 = np.array([10,20,30],int)
    need1 = np.array([15,20,25], int)
    new_matrix1 = NW.NW(matrix1,stock1,need1)
    stock1 = np.array([10,20,30],int)
    need1 = np.array([15,20,25], int)
    new_matrix2 = Min.Min(matrix1,stock1,need1)
    print(new_matrix1)
    print(new_matrix2)
