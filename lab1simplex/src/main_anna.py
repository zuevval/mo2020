import canon
import enter
import dual
from src import simplex, bruteforce
import numpy as np

if __name__ == "__main__":
    #original system
    A,b,c,func,MatrixSigns,VariablesSigns = enter.ReadFile('problem.txt')

    #canonical form
    CanonSys = canon.Convert(A,b,c,func,MatrixSigns, VariablesSigns)

    NpCanonForm = simplex.NpCanonicalForm(CanonSys)
    x = simplex.starting_vector(NpCanonForm)
    # x = simplex.starting_vector_method2(NpCanonForm).x  # TODO разобраться
    print("--- simplex algorithm: primal problem ---")
    print(simplex.simplex_method(NpCanonForm, x).x)
    print("------ brute force: primal problem ------")
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
    x = simplex.starting_vector(NpCanonForm)
    print("--- simplex algorithm: dual problem ---")
    print(simplex.simplex_method(NpCanonForm, x).x)
    print("------ brute force: dual problem ------")
    print(bruteforce(NpCanonForm))

