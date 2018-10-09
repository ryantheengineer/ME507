'''Code1.py: A Python script to complete coding assignment 1 for ME 507.'''

import numpy as np


# Establish a list of node numbers:
# n = [10, 100, 1000, 10000]
n = [10,100,1000]
# n = [10]

# Create structure of cases or iterations that go through the different cases
# where the definition of f changes

for elements in n:
    print("\n\nn = %d") % elements
    he = 1.0/elements
    # print("he = %f") % he

    # element-wise stiffness matrix
    ke = (1/he)*np.array([[1, -1],[-1, 1]])

    # print("ke = ", ke)

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
    F = np.zeros([elements,1])
    x1 = 0
    x2 = 0
    # Solving for the f elementwise vector
    # if f(x) = c:
    c = 1   # stand-in value for constant c
    for el in range(1,elements): # NOTE:May need to double check this range
        # print("el = %d") % el
        Ftemp = np.zeros([elements,1])  # reset Ftemp back to zeros
        x1 = (el-1)*he
        # print("x1 = %f") % x1
        x2 = el*he
        # print("x2 = %f") % x2
        fe[0] = (c/he)*((x2**2)/2 + x1*x2 + (x1**2)/2) # NOTE: consider making this a method to be called
        fe[1] = (c/he)*((x1**2)/2 - x1*x2 + (x2**2)/2)
        # print("fe = ",fe)

        Ftemp[el-1] = fe[0]
        Ftemp[el] = fe[1]
        # print("Ftemp = ",Ftemp)
        F += Ftemp
        # print("F = ",F)
    # NOTE: need to check to make sure the proper elements are being added into
    # the F vector at the end
    x1 = elements*he
    x2 = (elements+1)*he
    fe[0] = (c/he)*((x2**2)/2 + x1*x2 + (x1**2)/2)
    F[-1] += fe[0]
    # print("\n\n")
    # print("el = 10")
    # print("Adding stuff from element 10")
    print("F = ",F)

    d = np.linalg.solve(K,F)
    print("\n\n")
    print("d = ",d)
    # NOTE: this might be correct for the f(x) = c case
