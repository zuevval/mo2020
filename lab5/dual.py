import canon

def TransposeMatrix(A):
    B = []
    for _ in range(len(A[0])):
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

def  CreateDual(A,b,c,func,MatrixSigns,VariablesSigns):
    #dual system
    TransA = TransposeMatrix(A) #transposing matrix A
    Newb = NewVector(c) #finding dual vector b
    Newc = NewCoeff(b) #finding coefficients of dual goal function
    NewMatrixSigns = NewMSigns(VariablesSigns) #finding signs of dual matrix
    NewVariableSigns = NewVSigns(MatrixSigns) #finding signs of dual variables
    Newfunc = NewFunction(func) #finding dual problem type

    #canon dual system
    return canon.Convert(TransA,Newb,Newc,Newfunc,NewMatrixSigns,NewVariableSigns)


