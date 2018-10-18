'''Code1.py: A Python script to complete coding assignment 1 for ME 507.'''

import numpy as np
import matplotlib.pyplot as plt

# Linear shape function definitions
def N1(x1,x2,x):
    u = (x2 - x)/(x2 - x1)
    return u

def N2(x1,x2,x):
    u = (x - x1)/(x2 - x1)
    return u

def LM(a,e,n):        # a needs to be a 2x1 vector containing 1,2
    if a == 1:
        A = e
    if a == 2:
        A = e + 1
    if e == n and a == 2:   # NOTE: check this to make sure indices are handled correctly
        A = 0
    return A

# Establish a list of node numbers:
# n = [10, 100, 1000, 10000]
n = [10,100,1000]
# n = [10]


# Create structure of cases or iterations that go through the different cases
# where the definition of f changes
fcase = ['A','B','C']

c = 1

# Step through the different f(x) definitions
for i in range(0,3):    # NOTE: SHOULD BE (0,3)
    # for each # of elements in n vector
    print("\n")
    print("\n")
    print("fcase = ", fcase[i])

    # Step through the increasing numbers of elements for the chosen f(x)
    for elements in n:
        print("\n")
        print("n = ", elements)
        print("\n")

        N = 10*elements + 1
        x = np.linspace(0,1,N,endpoint=True)

        uA = np.zeros([len(x),1])
        uB = np.zeros([len(x),1])
        uC = np.zeros([len(x),1])
        for k in range(len(uA)):
            uA[k] = 0.5*c**2 - 0.5*c*x[k]**2
        for k in range(len(uB)):
            uB[k] = (1/6.)*(1 - x[k]**3)
        for k in range(len(uC)):
            uC[k] = (1/12.)*(1 - x[k]**4)

        u = np.zeros([N,3])
        for k in range(N):
            u[k,0] = uA[k]
            u[k,1] = uB[k]
            u[k,2] = uC[k]


        he = 1/float(elements)
        print("he = ",he)

        # Pull back Ke
        fe = np.zeros([2,1])
        K = np.zeros([elements,elements])
        F = np.zeros([elements,1])
        for e in range(1,elements + 1):
            Ke = np.zeros([2,2])
            for a in range(1,3):
                for b in range(1,3):
                    Ke[a-1][b-1] = (1/he)*(-1)**(a+b)
                    # Should result in indices from 0 to 1

            # Pull back fe
            f1 = 0
            f2 = 0
            c = 1
            x1 = he*(e-1)
            x2 = he*e
            if i == 0:
                f1 = c
                f2 = c
            elif i == 1:
                f1 = x1
                f2 = x2
            elif i == 2:
                f1 = x1**2
                f2 = x2**2
            fe[0] = (he/6)*(2*f1 + f2)
            fe[1] = (he/6)*(f1 + 2*f2)
            # print("fe = ",fe)

            # Assemble Ke and fe into K and F
            for a in range(1,3):
                if LM(a,e,elements) == 0:
                    continue
                else:
                    for b in range(1,3):
                        if LM(b,e,elements) == 0:
                            continue
                        else:
                            K[LM(a,e,elements)-1][LM(b,e,elements)-1] += Ke[a-1][b-1]

                F[LM(a,e,elements)-1] += fe[a-1]

        # print("F = ",F)
        # print("K = ",K)
        d = np.zeros([elements,1])
        d = np.linalg.solve(K,F)
        # np.append(d, [[0]], axis=0)
        # print("d = ",d)


        # create uh(x)
        uh = np.zeros(len(x))
        # print("uh = ",np.shape(uh))

        for e in range(1,elements+1):
            x1 = he*(e-1)
            x2 = he*e

            # set d1 and d2 for the given element
            d1 = d[e-1]
            if d1 == d[-1]:
                d2 = 0
            else:
                d2 = d[e]

            for j in range(len(x)):
                if x[j] >= x1 and x[j] < x2:
                    # print(x[j])
                    uh[j] += d1*N1(x1,x2,x[j]) + d2*N2(x1,x2,x[j])
                    # print(uh[j][i])
                else:
                    uh[j] += 0


        error = np.zeros(len(x))
        for m in range(len(error)):
            error[m] = uh[m]-u[m,i]
        # print("error = ",error)

        # Plot u(x) and uh(x) for the given combination of elements and load case
        title = ("FEA solution for load case " + fcase[i] + " with "
            + str(elements) + " elements")

        plt.figure(1)
        plt.title(title)
        plt.xlabel("Linear position along beam (x)")
        plt.ylabel("Displacement (u)")
        plt.plot(x,u[:,i],label='u(x)',linewidth=1,color='r')
        plt.plot(x,uh,label='uh(x)',linewidth=1,color='b',linestyle='-')
        # plt.plot(x,slope)
        plt.legend()
        plt.figure(2)
        plt.plot(x,error)
        # plt.show()

        ## Part 2: Compute global error ##
        # Given constants for computing error integral with 3-point Gauss quad.
        ksi = np.array([-np.sqrt(3./5.), 0, np.sqrt(3./5.)])
        w = np.array([5./9., 8./9., 5./9.])
        globerr = 0         # integral of error over entire beam
        dxksidksi = he/2    # constant required for integration
        uinterp = 0
        uhinterp = 0
        xe = 0
        d1 = 0
        d2 = 0

        for el in range(1,elements+1):
            for p in range(0,3):
                # We need u(xe(ksi[i])) and uh(ksi[i])
                x1 = he*(el-1)
                x2 = he*el
                # xpoints = [x1,x2]
                ksipoints = [-1,1]
                xe = (x1*N1(ksipoints[0],ksipoints[1],ksi[p])
                    + x2*N2(ksipoints[0],ksipoints[1],ksi[p]))
                # print("xe = ", xe)
                uinterp = np.interp(xe, x, u[:,p])

                d1 = d[el-1]
                if d1 == d[-1]:
                    d2 = 0
                else:
                    d2 = d[el]

                uhinterp = (d1*N1(ksipoints[0],ksipoints[1],ksi[p])
                        + d2*N2(ksipoints[0],ksipoints[1],ksi[p]))


                globerr += dxksidksi*w[p]*(np.abs(uinterp-uhinterp))**2
                # print("globerr = %f") % globerr
            # print("uinterp = ",uinterp)
            # print("uhinterp = ",uhinterp)
            # print("uinterp - uhinterp = ", uinterp-uhinterp)

        # print("globerr = %f") % globerr
        globerr = np.sqrt(globerr)
        print("globerr = %f") % globerr




        #
        # # Assemble the k element wise stiffness matrices into a global K matrix
        # K = np.zeros([elements,elements])
        # Ktemp = np.zeros([elements,elements])
        # # Here we need to map the parts of Ke onto the global matrix K
        # for j in range(0,elements-1):
        #     # Locate the Ke matrix in the global matrix
        #     Ktemp = np.zeros([elements,elements])
        #     Ktemp[j][j] = Ke[0][0]
        #     Ktemp[j][j+1] = Ke[0][1]
        #     Ktemp[j+1][j] = Ke[1][0]
        #     Ktemp[j+1][j+1] = Ke[1][1]
        #     # print("Ktemp = ",Ktemp)
        #
        #     # Add Ktemp into the global matrix
        #     K += Ktemp
        #     # print("K = ", K)
        # # Add the final piece to the bottom right element of K:
        # K[elements-1][elements-1] += Ke[0][0]
        # # print("K = ",K)
        #
        # # Find F
        # fe = np.zeros([2,1])
        # F = np.zeros([elements,1])
        # x1 = 0
        # x2 = 0
        #
        # # cycle through elements and solve for fe
        # for el in range(1,elements+1):
        #     Ftemp = np.zeros([elements,1])
        #
        #     # print("element # %d") % el
        #     x1 = he*(el-1)
        #     x2 = he*el
        #
        #     if i == 0:
        #         fe[0] = float(c*he/2)
        #         fe[1] = float(c*he/2)
        #
        #     if i == 1:
        #         fe[0] = (1/he)*((x1**3)/3 - (x2*x1**2)/2 + (x2**3)/6)
        #         fe[1] = (1/he)*((x1**3)/6 - (x1*x2**2)/2 + (x2**3)/3)
        #
        #     if i == 2:
        #         fe[0] = (1/he)*((x1**4)/4 - (x2*x1**3)/3 + (x2**4)/12)
        #         fe[1] = (1/he)*((x1**4)/12 - (x1*x2**3)/3 + (x2**4)/4)
        #
        #
        #     if el == elements:
        #         Ftemp[-1] = fe[0]
        #     else:
        #         Ftemp[el-1] = fe[0]
        #         Ftemp[el] = fe[1]
        #
        #     # print("Ftemp = ",Ftemp)
        #     F += Ftemp
        #
        # # print("F = ",F)
        #
        # # Find d
        # d = np.zeros([elements,1])
        # d = np.linalg.solve(K,F)
        # print("d = ",d)
        #
        #
        # # create uh(x)
        # uh = np.zeros(len(x))
        # # print("uh = ",np.shape(uh))
        #
        # for el in range(1,elements+1):
        #     x1 = he*(el-1)
        #     x2 = he*el
        #
        #     # set d1 and d2 for the given element
        #     d1 = d[el-1]
        #     if d1 == d[-1]:
        #         d2 = 0
        #     else:
        #         d2 = d[el]
        #
        #     for j in range(len(x)):
        #         if x[j] >= x1 and x[j] < x2:
        #             # print(x[j])
        #             uh[j] += d1*N1(x1,x2,x[j]) + d2*N2(x1,x2,x[j])
        #             # print(uh[j][i])
        #         else:
        #             uh[j] += 0
        #
        # # slope = np.zeros(len(uh))
        # # for m in range(len(slope)-1):
        # #     slope[m] = (uh[m+1]-uh[m])/(x[m+1]-x[m])
        # # print("uh = ",np.shape(uh))
        # # print(slope)
        # error = np.zeros(len(x))
        # for m in range(len(error)):
        #     error[m] = uh[m]-u[m,i]
        # # print("error = ",error)
        #
        # # Plot u(x) and uh(x) for the given combination of elements and load case
        # title = ("FEA solution for load case " + fcase[i] + " with "
        #     + str(elements) + " elements")
        #
        # plt.figure(1)
        # plt.title(title)
        # plt.xlabel("Linear position along beam (x)")
        # plt.ylabel("Displacement (u)")
        # plt.plot(x,u[:,i],label='u(x)',linewidth=1,color='r')
        # plt.plot(x,uh,label='uh(x)',linewidth=1,color='b',linestyle='-')
        # # plt.plot(x,slope)
        # plt.legend()
        # # plt.figure(2)
        # # plt.plot(x,error)
        # # plt.show()
        #
        #
        # ## Part 2: Compute global error ##
        # # Given constants for computing error integral with 3-point Gauss quad.
        # ksi = np.array([-np.sqrt(3./5.), 0, np.sqrt(3./5.)])
        # w = np.array([5./9., 8./9., 5./9.])
        # globerr = 0         # integral of error over entire beam
        # dxksidksi = he/2    # constant required for integration
        # uinterp = 0
        # uhinterp = 0
        # xe = 0
        # d1 = 0
        # d2 = 0
        #
        # for el in range(1,elements+1):
        #     for p in range(0,3):
        #         # We need u(xe(ksi[i])) and uh(ksi[i])
        #         x1 = he*(el-1)
        #         x2 = he*el
        #         # xpoints = [x1,x2]
        #         ksipoints = [-1,1]
        #         xe = (x1*N1(ksipoints[0],ksipoints[1],ksi[p])
        #             + x2*N2(ksipoints[0],ksipoints[1],ksi[p]))
        #         # print("xe = ", xe)
        #         uinterp = np.interp(xe, x, u[:,p])
        #
        #         d1 = d[el-1]
        #         if d1 == d[-1]:
        #             d2 = 0
        #         else:
        #             d2 = d[el]
        #
        #         uhinterp = (d1*N1(ksipoints[0],ksipoints[1],ksi[p])
        #                 + d2*N2(ksipoints[0],ksipoints[1],ksi[p]))
        #
        #
        #         globerr += dxksidksi*w[p]*(np.abs(uinterp-uhinterp))**2
        #         # print("globerr = %f") % globerr
        #     # print("uinterp = ",uinterp)
        #     # print("uhinterp = ",uhinterp)
        #     # print("uinterp - uhinterp = ", uinterp-uhinterp)
        #
        # # print("globerr = %f") % globerr
        # globerr = np.sqrt(globerr)
        # print("globerr = %f") % globerr
        #
        # ##
        # # plt.figure(2)
        # # plt.loglog(e,e)
