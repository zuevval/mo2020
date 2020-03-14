
def NW(matrix, stock, need):
    new_matrix = []
    for _ in range(len(matrix)):
        l = []
        for _ in range(len(matrix[0])):
            l.append(-1)                
            #-1 defines empty cell
        new_matrix.append(l)
    
    i = 0
    j = 0
    while i != len(matrix) and j != len(matrix[0]):
        if stock[i] < need[j]:
            new_matrix[i][j] = stock[i]
            need[j] = need[j] - stock[i]
            stock[i] = 0
            i = i + 1
        elif stock[i] > need[j]:
            new_matrix[i][j] = need[j]
            stock[i] = stock[i] - need[j]
            need[j] = 0
            j = j + 1
        else:
            new_matrix[i][j] = need[j]
            stock[i] = 0
            need[j] = 0
            i = i + 1
            j = j + 1
    return new_matrix


