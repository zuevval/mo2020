import numpy as np
MAX = pow(10,6)

def Min(matrix, stock, need):
    new_matrix = []
    help_matrix = matrix.view()
    for _ in range(len(matrix)):
        l = []
        for _ in range(len(matrix[0])):
            l.append(-1)                
            #-1 defines empty cell
        new_matrix.append(l)
    
    counter = 0
    while len(matrix) + len(matrix[0]) - 1 != counter:
        min_str = list(min(i) for i in help_matrix)
        min_elements = np.where(help_matrix ==  min(min_str))
        
        i,j = min_elements[0][0], min_elements[1][0]

        help_matrix[i][j] = MAX
        if stock[i] == 0 or need[j] == 0:
            continue
        counter = counter + 1

        if stock[i] < need[j]:
            new_matrix[i][j] = stock[i]
            need[j] = need[j] - stock[i]
            stock[i] = 0
        elif stock[i] > need[j]:
            new_matrix[i][j] = need[j]
            stock[i] = stock[i] - need[j]
            need[j] = 0
        else:
            new_matrix[i][j] = need[j]
            stock[i] = 0
            need[j] = 0

    return new_matrix


