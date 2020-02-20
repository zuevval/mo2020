import canon
import dual
from src import simplex
from bruteforce import bruteforce

if __name__ == "__main__":
    #original system
    A = [[-4,2,3], [-2,-2,1], [-1,-1,3]] 
    b = [8,2,3]
    c = [2,0,-1] #coefficiens for goal function
    func = "max" #type of a problem
    MatrixSigns = ["=", "<=", ">="] #matrix restrictions signs
    VariablesSigns = ["<=", ">=", ""] # "" if no restictions on x

    #canonical form
    CanonSys = canon.Convert(A,b,c,func,MatrixSigns, VariablesSigns)

    NpCanonForm = simplex.NpCanonicalForm(CanonSys)
    x = simplex.starting_vector(NpCanonForm)
    stopIteration = False
    while not stopIteration:
        x, stopIteration = simplex.simplex_step(NpCanonForm, x)
        print(simplex.np.dot(x, NpCanonForm.c))

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
