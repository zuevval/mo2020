import canon
import dual
from src import simplex, bruteforce
import numpy as np

if __name__ == "__main__":
    #original system
    A = [[1,-2,2,0], [1,2,1,1], [2,1,-4,1], [-1,4,0,-2]] 
    b = [6,24,30,-6]
    c = [8,3,4,2] #coefficiens for goal function
    func = "max" #type of a problem
    MatrixSigns = ["<=", "=", "=", ">="] #matrix restrictions signs
    VariablesSigns = [">=", "", "", ""] # "" if no restictions on x

    #canonical form
    CanonSys = canon.Convert(A,b,c,func,MatrixSigns, VariablesSigns)

    NpCanonForm = simplex.NpCanonicalForm(CanonSys)
    x = simplex.starting_vector(NpCanonForm)
    # x = simplex.starting_vector_method2(NpCanonForm).x  # TODO разобраться
    print(simplex.simplex_method(NpCanonForm, x).x)

    print(bruteforce(NpCanonForm))

    #dual system
    dual.ConvertToMax(c, func) #always solving ->max
    dual.CheckSigns(A, b, MatrixSigns) #only <= sign allowed in finding max problem
    TransA = dual.TransposeMatrix(A) #transposing matrix A
    Newb = dual.NewVector(c) #finding dual vector b
    Newc = dual.NewCoeff(b) #finding coefficients of dual goal function
    NewMatrixSigns = dual.NewMSigns(VariablesSigns) #finding signs of dual matrix
    NewVariableSigns = dual.NewVSigns(MatrixSigns) #finding signs of dual variables
    Newfunc = dual.NewFunction(func) #finding dual problem type

    #canon dual system
    DualCanonSys = canon.Convert(TransA,Newb,Newc,func,NewMatrixSigns,NewVariableSigns)

    NpCanonForm = simplex.NpCanonicalForm(DualCanonSys)
    print("A: "+ str(NpCanonForm.A))
    print("rank(A): " + str(np.rank(NpCanonForm.A)))
    x = simplex.starting_vector(NpCanonForm)
    print(simplex.simplex_method(NpCanonForm, x).x)

    print(bruteforce(NpCanonForm))

