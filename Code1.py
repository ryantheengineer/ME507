'''Code1.py: A Python script to complete coding assignment 1 for ME 507.'''

import numpy as np
import matplotlib.pyplot as plt
import pyprind

# Linear shape function definitions
def N1(x1,x2,x):
    u = (x2 - x)/(x2 - x1)
    return u

def N2(x1,x2,x):
    u = (x - x1)/(x2 - x1)
    return u

# Establish a list of node numbers:
# n = [10, 100, 1000, 10000]
n = [10,100,1000]
# n = [10]


# Create structure of cases or iterations that go through the different cases
# where the definition of f changes
fcase = ['A','B','C']
c = 1

N = 100
x = np.linspace(0,1,N,endpoint=True)
uA = np.zeros([len(x),1])
uB = np.zeros([len(x),1])
uC = np.zeros([len(x),1])
for i in range(len(uA)):
    uA[i] = 0.5*c**2 - 0.5*c*x[i]**2
for i in range(len(uB)):
    uB[i] = (1/6.)*(1 - x[i]**3)
for i in range(len(uC)):
    uC[i] = (1/12.)*(1 - x[i]**4)

u = np.zeros([N,3])
for i in range(N):
    u[i,0] = uA[i]
    u[i,1] = uB[i]
    u[i,2] = uC[i]



for i in range(0,2):    # NOTE: SHOULD BE (0,3)
    # for each # of elements in n vector
    print("\n")
    print("fcase = ", fcase[i])
    print("\n")

    for elements in n:
        print("\n")
        print("elements = ", elements)
        print("\n")

        # Find K
        he = 1/float(elements)
        ke = (1/he)*np.array([[1, -1],[-1, 1]])
        # Assemble the k element wise stiffness matrices into a global K matrix
        K = np.zeros([elements,elements])
        Ktemp = np.zeros([elements,elements])
        # Here we need to map the parts of ke onto the global matrix K
        for j in range(0,elements-1):
            # Locate the ke matrix in the global matrix
            Ktemp = np.zeros([elements,elements])
            Ktemp[j][j] = ke[0][0]
            Ktemp[j][j+1] = ke[0][1]
            Ktemp[j+1][j] = ke[1][0]
            Ktemp[j+1][j+1] = ke[1][1]
            # print("Ktemp = ",Ktemp)

            # Add Ktemp into the global matrix
            K += Ktemp
            # print("K = ", K)
        # Add the final piece to the bottom right element of K:
        K[elements-1][elements-1] += ke[0][0]
        # print("K = ",K)

        # Find F
        fe = np.zeros([2,1])
        F = np.zeros([elements,1])
        x1 = 0
        x2 = 0

        # cycle through elements and solve for fe
        for el in range(1,elements+1):
            Ftemp = np.zeros([elements,1])

            # print("element # %d") % el
            x1 = he*(el-1)
            x2 = he*el

            if i == 0:
                fe[0] = float(c*he/2)
                fe[1] = float(c*he/2)

            if i == 1:
                fe[0] = (1/he)*((x1**3)/3 - (x2*x1**2)/2 + (x2**3)/6)
                fe[1] = (1/he)*((x1**3)/6 - (x1*x2**2)/2 + (x2**3)/3)

            if i == 2:
                fe[0] = 2
                fe[1] = 2


            if el == elements:
                Ftemp[-1] = fe[0]
            else:
                Ftemp[el-1] = fe[0]
                Ftemp[el] = fe[1]

            # print("Ftemp = ",Ftemp)
            F += Ftemp

        print("F = ",F)

        # Find d
        d = np.zeros([elements,1])
        d = np.linalg.solve(K,F)
        print("d = ",d)


        # create uh(x)
        # graph u(x) and uh(x) for the given combination of elements and load case






















#     print("\nf(x) for case %s") % loadcase
#     for elements in n:
#         print("\nn = %d") % elements
#         he = 1.0/elements
#
#         # element-wise stiffness matrix
#         ke = (1/he)*np.array([[1, -1],[-1, 1]])
#
#         # Assemble the k element wise stiffness matrices into a global K matrix
#         K = np.zeros([elements,elements])
#         Ktemp = np.zeros([elements,elements])
#         # Here we need to map the parts of ke onto the global matrix K
#         for i in range(0,elements-1):
#             # Locate the ke matrix in the global matrix
#             Ktemp = np.zeros([elements,elements])
#             Ktemp[i][i] = ke[0][0]
#             Ktemp[i][i+1] = ke[0][1]
#             Ktemp[i+1][i] = ke[1][0]
#             Ktemp[i+1][i+1] = ke[1][1]
#             # print("Ktemp = ",Ktemp)
#
#             # Add Ktemp into the global matrix
#             K += Ktemp
#             # print("K = ", K)
#         # Add the final piece to the bottom right element of K:
#         K[elements-1][elements-1] += ke[0][0]
#         print("K = ",K)
#         if loadcase == 'A' and elements == 10:
#             KA10 = K
#         if loadcase == 'B' and elements == 10:
#             KB10 = K
#         if loadcase == 'C' and elements == 10:
#             KC10 = K
#         if loadcase == 'A' and elements == 100:
#             KA100 = K
#         if loadcase == 'B' and elements == 100:
#             KB100 = K
#         if loadcase == 'C' and elements == 100:
#             KC100 = K           # NOTE: add code for 1000 and up elements later
#
#         # NOTE: The n = 10000 node problem may be solved by simplifying K down in such
#         # a way that only the nonzero indices are used.
#
#         fe = np.zeros([2,1])
#         F = np.zeros([elements,1])
#         x1 = 0
#         x2 = 0
#         # Solving for the f elementwise vector
#         # if f(x) = c:
#         c = 1  # stand-in value for constant c
#         for el in range(1,elements): # NOTE:May need to double check this range
#             # print("el = %d") % el
#             Ftemp = np.zeros([elements,1])  # reset Ftemp back to zeros
#             x1 = (el-1)*he
#             # print("x1 = %f") % x1
#             x2 = el*he
#             # print("x2 = %f") % x2
#             if loadcase == 'A':
#                 fe[0] = (c/he)*((x2**2)/2 + x1*x2 + (x1**2)/2) # NOTE: consider making this a method to be called
#                 fe[1] = (c/he)*((x1**2)/2 - x1*x2 + (x2**2)/2)
#             elif loadcase == 'B':
#                 fe[0] = (1/he)*((x1**3)/3 - (x2*x1**2)/2 + (x2**3)/6)
#                 fe[1] = (1/he)*((x1**3)/6 - (x1*x2**2)/2 + (x2**3)/3)
#             elif loadcase == 'C':
#                 fe[0] = (1/he)*((x1**4)/4 - (x2*x1**3)/3 + (x2**4)/12)
#                 fe[1] = (1/he)*((x1**4)/12 - (x1*x2**3)/3 + (x2**4)/4)
#             # print("fe = ",fe)
#
#             Ftemp[el-1] = fe[0]
#             Ftemp[el] = fe[1]
#             # print("Ftemp = ",Ftemp)
#             F += Ftemp
#             # print("F = ",F)
#         # NOTE: need to check to make sure the proper elements are being added into
#         # the F vector at the end
#         x1 = elements*he
#         x2 = (elements+1)*he
#         fe[0] = (c/he)*((x2**2)/2 + x1*x2 + (x1**2)/2)
#         F[-1] += fe[0]
#         # print("\n\n")
#         # print("el = 10")
#         # print("Adding stuff from element 10")
#         print("F = ",F)
#         if loadcase == 'A':
#             FA = F
#         if loadcase == 'B':
#             FB = F
#         if loadcase == 'C':
#             FC = F
#
#         d = np.linalg.solve(K,F)
#         if loadcase == 'A':
#             dA = d
#         if loadcase == 'B':
#             dB = d
#         if loadcase == 'C':
#             dC = d
#         # print(np.shape(d))
#         print("d = ",d)
#         # NOTE: this might be correct for the f(x) = c case
#
# # PLOT THE EXACT AND APPROXIMATE SOLUTION FOR EACH F(X) AND N:
#
# # Create xvector for exact solutions (smooth curve)
# N = 100
# x = np.linspace(0,1,N,endpoint=True)
# uA = np.zeros([len(x),1])
# uB = np.zeros([len(x),1])
# uC = np.zeros([len(x),1])
# for i in range(len(uA)):
#     uA[i] = 0.5*c**2 - 0.5*c*x[i]**2
# for i in range(len(uB)):
#     uB[i] = (1/6.)*(1 - x[i]**3)
# for i in range(len(uC)):
#     uC[i] = (1/12.)*(1 - x[i]**4)
#
# u = np.zeros([N,3])
# for i in range(N):
#     u[i,0] = uA[i]
#     u[i,1] = uB[i]
#     u[i,2] = uC[i]
# # print("uA = ",uA)
# # print("uB = ",uB)
# print("u = ",u)
# print("\n")
# print("size of u is ",np.shape(u))
# print("\n")
#
# # u =  [uA,uB,uC]
# uh = np.zeros([len(x),3])
# # he = [heA,heB,heC]
# # print("he = ", he)
#
# # NOTE: HAVE K MATRICES BE CALLED SPECIFICALLY, BASED ON WHICH MATRIX AND N COMBO
# # IS APPROPRIATE FOR THE GIVEN CALCULATION?
#
#
# # K10 = np.zeros([10,10,3])
# # for i in range(10):
# #     for j in range(10):
# #         K10[i][j][0] = KA10[i][j]
# #         K10[i][j][1] = KB10[i][j]
# #         K10[i][j][2] = KC10[i][j]
# # # K10 =  [KA10,KB10,KC10]
# # print("K10 = ", K10)
# # print("KA10 = ", K10[:][:][0])
#
# F =  [FA,FB,FC]
# d = np.zeros([10,3])
# for i in range(10):
#     d[i,0] = dA[i]
#     d[i,1] = dB[i]
#     d[i,2] = dC[i]
#
# # d =  [dA,dB,dC]
# # print(np.shape(d))
# # print(d)
#
# # m = len(x)
# # bar = pyprind.ProgBar(m)
# for i in range(len(fcase)):     # NOTE: This may be the wrong range to iterate over
#     # print(range(len(fcase)))
#     # PSEUDO CODE FOR APPROXIMATE SOLUTIONS:
#     # set the number of elements
#     # step through the elements
#     for elements in n:
#         he = 1.0/elements
#         for el in range(1,elements+1):    # NOTE: make sure to account for d11 (end point of element 10) being zero (due to constraint)
#             # for each element, set the endpoint x1 and x2 values
#             x1 = (el-1)*he
#             x2 = el*he
#             # print("x1 = %f") % x1
#             # print("x2 = %f") % x2
#             # set d1 and d2 for the given element
#             d1 = d[el-1][i]
#             if d1 == d[-1][i]:
#                 d2 = 0
#             else:
#                 d2 = d[el][i]
#             # print("d1 = %f") %d1
#             # print("d2 = %f") %d2
#
#             # calculate uh for the given n over the applicable x range
#             for j in range(len(x)):
#                 if x[j] >= x1 and x[j] < x2:
#                     # print(x[j])
#                     uh[j][i] += d1*N1(x1,x2,x[j]) + d2*N2(x1,x2,x[j])
#                     # print(uh[j][i])
#                 else:
#                     uh[j][i] += 0
#                     # print(uh[j][i])
#                 # bar.update()
# # print("uh = ", uh)
#         x1 = (elements-1)*he
#         x2 = elements*he
#         # print("x1 = %f") % x1
#         # print("x2 = %f") % x2
#
#         d1 = d[elements-1][i]
#         d2 = 0
#         # print("d1 = %f") %d1
#         # print("d2 = %f") %d2
#
#         for j in range(len(x)):
#             if x[j] >= x1 and x[j] < x2:
#                 # print(x[j])
#                 uh[j][i] += d1*N1(x1,x2,x[j]) + d2*N2(x1,x2,x[j])
#                 # print(uh[j][i])
#             else:
#                 uh[j][i] += 0
#
# print("uh = ", uh)
# print("\n")
# print("size of uh is ",np.shape(uh))
# print("\n")
#
#
#     # assemble the uh values for the element into a larger uh vector, ready to plot against x later
#
#
# plt.figure()
# plt.plot(x,u[:,0])
# # plt.plot(x,uh[:,0])
# # plt.plot(x,uh[:,1])
# # plt.plot(x,uh[:,2])
# plt.show()
#
#     # for elements in n:
#     #     he = 1.0/elements
#     #     for el in range(1,elements+1):
#     #         # set the element number el
#     #         # determine endpoints of the given element
#     #         x1 = (el-1)*he
#     #         x2 = el*he
#     #         for dnum in d:
#     #             uh = dnum*N1(x1,x2,x)
#     #         while x >= x1 and x < x2:
#     #             plt.plot(x,)
