def ReadFile(file_name):
    f = open(file_name)
    A = []
    for line in f:
        if line == "\n":
            break
        A.append([int(x) for x in line.split()])

    for line in f:
        if line == "\n":
            break
        MatrixSigns = [x for x in line.split()]
    
    for line in f:
        if line == "\n":
            break
        b = [int(x) for x in line.split()]

    for line in f:
        if line == "\n":
            break
        c = [int(x) for x in line.split()]


    for line in f:
        if line == "\n":
            break
        VariablesSigns = [x for x in line.split()]

    for _ in range(len(c)):
                VariablesSigns.append("")

    for line in f:
        if line == "\n":
            break
        func = line

    print("A = ", A)
    print("MatrixSigns = ", MatrixSigns)
    print("b = ", b)
    print("c = ", c)
    print("VariablesSigns = ", VariablesSigns)
    print("func = ", func)
    return A,b,c,func,MatrixSigns,VariablesSigns

