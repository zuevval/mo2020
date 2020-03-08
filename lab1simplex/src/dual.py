def CheckSigns(A,b,MatrixSigns):
    for i in range(len(MatrixSigns)):
        if MatrixSigns[i] == "<=":
            for j in range(len(A[i])):
                A[i][j]=-A[i][j]
            b[i]= -b[i]
            MatrixSigns[i] = ">="

def TransposeMatrix(A):
    B = []
    for _ in range(len(A)):
        B.append([])
    for i in range(len(A)):
        for j in range (len(A[i])):
            B[j].append(A[i][j])
    return B

def NewVector(c):
    return c

def NewCoeff(b):
    return b

def NewMSigns(VariablesSigns):
    signs = []
    for i in range(len(VariablesSigns)):
        if VariablesSigns[i] == ">=":
            signs.append(">=")
        elif VariablesSigns[i] == "":
            signs.append("=")
    return signs

def NewVSigns(MatrixSigns):
    signs = []
    for i in range(len(MatrixSigns)):
        if MatrixSigns[i] == ">=":
            signs.append(">=")
        elif MatrixSigns[i] == "=":
            signs.append("")
    return signs

def NewFunction(func):
    if func == "min":
        return "max"
    return "min"



