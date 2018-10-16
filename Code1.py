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

# Step through the different f(x) definitions
for i in range(0,3):    # NOTE: SHOULD BE (0,3)


    # for each # of elements in n vector
    print("\n")
    print("\n")
    print("fcase = ", fcase[i])

    # Step through the increasing numbers of elements for the chosen f(x)
    for elements in n:
        print("\n")
        print("elements = ", elements)
        print("\n")

        N = 10*elements
        x = np.linspace(0,1,N,endpoint=False)

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
                fe[0] = (1/he)*((x1**4)/4 - (x2*x1**3)/3 + (x2**4)/12)
                fe[1] = (1/he)*((x1**4)/12 - (x1*x2**3)/3 + (x2**4)/4)


            if el == elements:
                Ftemp[-1] = fe[0]
            else:
                Ftemp[el-1] = fe[0]
                Ftemp[el] = fe[1]

            # print("Ftemp = ",Ftemp)
            F += Ftemp

        # print("F = ",F)

        # Find d
        d = np.zeros([elements,1])
        d = np.linalg.solve(K,F)
        # print("d = ",d)


        # create uh(x)
        uh = np.zeros(len(x))
        # print("uh = ",np.shape(uh))

        for el in range(1,elements+1):
            x1 = he*(el-1)
            x2 = he*el

            # set d1 and d2 for the given element
            d1 = d[el-1]
            if d1 == d[-1]:
                d2 = 0
            else:
                d2 = d[el]

            for j in range(len(x)):
                if x[j] >= x1 and x[j] < x2:
                    # print(x[j])
                    uh[j] += d1*N1(x1,x2,x[j]) + d2*N2(x1,x2,x[j])
                    # print(uh[j][i])
                else:
                    uh[j] += 0

        # for el in range(1,elements+1):
        #     x1 = he*el(-1)

        # slope = np.zeros(len(uh))
        # for m in range(len(slope)-1):
        #     slope[m] = (uh[m+1]-uh[m])/(x[m+1]-x[m])
        # print("uh = ",np.shape(uh))
        # print(slope)
        error = np.zeros(len(x))
        for m in range(len(error)):
            error[m] = uh[m]-u[m,i]

        # Plot u(x) and uh(x) for the given combination of elements and load case
        title = ("FEA solution for load case " + fcase[i] + " with "
            + str(elements) + " elements")

        # plt.figure(1)
        # plt.title(title)
        # plt.xlabel("Linear position along beam (x)")
        # plt.ylabel("Displacement (u)")
        # plt.plot(x,u[:,i],label='u(x)',linewidth=1,color='r')
        # plt.plot(x,uh,label='uh(x)',linewidth=1,color='b',linestyle='-')
        # # plt.plot(x,slope)
        # plt.legend()
        # plt.figure(2)
        # plt.plot(x,error)
        # plt.show()


        ## Part 2: Compute global error ##
        # Given constants for computing error integral with 3-point Gauss quad.
        ksi = np.array([-np.sqrt(3./5.), 0, np.sqrt(3./5.)])
        w = np.array([5./9., 8./9., 5./9.])
        globerr = 0         # integral of error over entire beam
        dxksidksi = he/2    # constant required for integration
        uinterp = 0
        uhinterp = 0

        for el in range(1,elements+1):
            for p in range(0,3):
                # Need to interpolate the value for u and uh at the current
                # value of ksi to give to the integration function.
                # print('shape of x: ',np.shape(x))
                # print('shape of u: ',np.shape(u))
                # print('shape of uh: ',np.shape(uh))
                x1 = he*(el-1)
                x2 = he*el
                # xpoints = [x1,x2]
                # ksipoints = [-1,1]
                # Given the ksi values (which don't change), map those to the
                # appropriate x value found in the given element
                ksiinterp = 0
                xinterp = 0
                ksiinterp = (ksi[p] + 1)/2. # Normalize ksi value on parent domain
                # ksi in geometric domain, ready for interpolation
                xinterp = x1 + ksiinterp*(x2 - x1)

                uinterp = np.interp(xinterp, x, u[:,p])
                uhinterp = np.interp(xinterp, x, uh)

                globerr += dxksidksi*w[p]*(np.abs(uinterp-uhinterp))**2
                # print("globerr = %f") % globerr
            # print("uinterp = ",uinterp)
            # print("uhinterp = ",uhinterp)
            # print("uinterp - uhinterp = ", uinterp-uhinterp)

        print("globerr = %f") % globerr
        globerr = np.sqrt(globerr)
        print("globerr = %f") % globerr
