import utils as cf
import enter
def ConvertToEquals(MatrixSigns):
    w = []
    for i in range(len(MatrixSigns)):
        if MatrixSigns[i] == "<=":
            w.append(1)
        elif MatrixSigns[i] == ">=":
            w.append(-1)
        else:
            w.append(0)
    return w

def CheckVariablesSigns(VariablesSigns, x):
    u = []
    for i in range(len(VariablesSigns)):
        if VariablesSigns[i] == "<=":
            x[i] = -x[i]
            VariablesSigns[i] = ">="
            u.append(0)
        elif VariablesSigns[i] == "":
            u.append(1)
        else:
            u.append(0)
    return u
    
def PrintSystem(A,b,c,x,vectors,w,func, MatrixSigns,VariablesSigns):
    for i in range(len(c)):
        if c[i] != 0:
            if vectors != [] and vectors[i] != 0:
                print(c[i], "(u[", i, "]", "-", "v[", i,"])", sep = '', end = '')
            else:
                print(c[i], "x[", i, "]", sep = '', end = '')
        if i != len(c) - 1 and c[i+1]>0:
            print("+", end = '')
    if w != []:
        for i in range(len(w)):
            if w[i] != 0:
                print("+0w[", i, "]", sep = '', end = '')
    print(" ->", func)
    for i in range(len(A)):
        for j in range(len(A[i])):   
            if vectors != [] and vectors[j] != 0:
                print(A[i][j], "(u[", j, "]", "-", "v[", j,"])", sep = '', end = '') 
            else:
                print(A[i][j], "x[", j, "]", sep = '', end = '') 
            if j != len(A[i]) - 1:
                if A[i][j + 1] >= 0:
                    print("+", sep = '', end = '')
        if w != []:
            if w[i] > 0:
                print("+", end = '')
            if w[i] != 0:
                print(w[i], "w[", i,"]", "=", b[i], sep = '')
            else:
                print(MatrixSigns[i], b[i], sep = '')
        else:
            print(MatrixSigns[i], b[i], sep = '')
    for i in range(len(VariablesSigns)):
        if VariablesSigns[i] != "":
            if x[i] > 0:
                print("x[", i, "]", VariablesSigns[i], "0", sep = '', end = ' ')
            else:
                print("-x[", i, "]", VariablesSigns[i], "0", sep = '', end = ' ')
    for i in range(len(vectors)):
        if vectors[i] != 0:
            print("u[", i, "]>=0", sep = '', end = ' ')
            print("v[", i, "]>=0", sep = '', end = ' ')     
    for i in range(len(w)):
        if w[i] != 0:
            print("w[", i, "]>=0", sep = '', end = ' ')
    print(end = '\n\n')

def CreateSystem(A,b,c,vectors,w):
    B = []
    k = 0

    for _ in range(len(A)):
        B.append([])

    j = 0
    for i in range(len(B)):
        for j in range(len(A[i])):
            B[i].append(A[i][j])
            if vectors != [] and vectors[j] != 0:
                B[i].append(-A[i][j])

    d = []
    for i in range(len(c)):
        d.append(c[i])
        if vectors != [] and vectors[i] != 0:
            d.append(-c[i])

    
    if w != []:
        for i in range(len(w)):
            if w[i] != 0:
                for j in range(len(B)):
                    if i == j:
                        B[j].append(w[i])
                    else:
                        B[j].append(0)
        for i in range(len(w)):
            if w[i] != 0:
                d.append(0)

    return cf.CanonicalForm(B,b,d)   
            
def Convert(A,b,c,func,MatrixSigns,VariablesSigns):
    x = []
    for i in range(len(c)):
        x.append(1)
    PrintSystem(A,b,c,x,[],[],func,MatrixSigns,VariablesSigns) #print given system
    w = ConvertToEquals(MatrixSigns) #adding vector w to get rid off <= and >= and make them =
    vectors = CheckVariablesSigns(VariablesSigns, x) #adding vectors u and v to set
    #restriction >=0 on variables that can be of any sign
    PrintSystem(A,b,c,x,vectors,w,func,MatrixSigns,VariablesSigns) #print converted to canon system
    return CreateSystem(A,b,c,vectors,w)


A,b,c,func,MatrixSigns,VariablesSigns = enter.ReadFile('problem.txt')
Convert(A,b,c,func,MatrixSigns, VariablesSigns)
