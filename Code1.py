'''Code1.py: A Python script to complete coding assignment 1 for ME 507.'''

import numpy as np


# Establish a list of node numbers:
# n = [10, 100, 1000, 10000]
# n = [10,100,1000]
n = [10]

# Create structure of cases or iterations that go through the different cases
# where the definition of f changes

for elements in n:
    he = 1.0/elements

    # element-wise stiffness matrix
    ke = (1/he)*np.array([[1, -1],[-1, 1]])

    print("ke = ", ke)

    # Assemble the k element wise stiffness matrices into a global K matrix
    K = np.zeros([elements,elements])
    Ktemp = np.zeros([elements,elements])
    # Here we need to map the parts of ke onto the global matrix K
    for i in range(0,elements-1):
        # Locate the ke matrix in the global matrix
        Ktemp = np.zeros([elements,elements])
        Ktemp[i][i] = ke[0][0]
        Ktemp[i][i+1] = ke[0][1]
        Ktemp[i+1][i] = ke[1][0]
        Ktemp[i+1][i+1] = ke[1][1]
        # print("Ktemp = ",Ktemp)

        # Add Ktemp into the global matrix
        K += Ktemp
        # print("K = ", K)
    # Add the final piece to the bottom right element of K:
    K[elements-1][elements-1] += ke[0][0]
    print("K = ",K)

    # NOTE: The n = 10000 node problem may be solved by simplifying K down in such
    # a way that only the nonzero indices are used.

    fe = np.zeros([2,1])
    x1 = 0
    x2 = 0
    # Solving for the f elementwise vector
    # if f(x) = c:
    c = 1   # stand-in value for constant c
    for el in range(1,elements+1): # NOTE:May need to double check this range
        x1 = (el-1)*he
        x2 = el*he
        fe[0] = (c/he)*((x2**2)/2 + x1*x2 + (x1**2)/2) # NOTE: consider making this a method to be called
        fe[1] = (c/he)*((x1**2)/2 - x1*x2 + (x2**2)/2)
        print("fe = ",fe)

        # NOTE: Then need to add current fe vector into the proper place in Ftemp,
        # which then should be added into the global F vector, like the method
        # for K
        
