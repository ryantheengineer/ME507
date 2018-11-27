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
    if e == n and a == 2:
        A = 0
    return A

# Establish a list of node numbers:
# n = [10, 100, 1000, 10000]
n = [10,100,1000]
# n = [10]

# Create structure of cases or iterations that go through the different cases
# where the definition of f changes
fcase = ['A','B','C']

c = 1   # Arbitrary constant for load case A

# Step through the different f(x) definitions
for i in range(0,3):
    print('\n')
    print('\n')
    print('fcase = ', fcase[i])
    evector = np.zeros([len(n),1])
    count = 0
    # Step through the increasing numbers of elements for the chosen f(x)
    for elements in n:
        print('\n')
        print('n = ', elements)
        print('\n')
        N = 10
        uh = np.zeros([1,1])
        uhtemp = np.zeros([N,1])
        x = np.zeros([1,1])
        xtemp = np.zeros([N,1])
        he = 1/float(elements)

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

        d = np.zeros([elements,1])
        d = np.linalg.solve(K,F)
        # print('K = ',K)
        # print('F = ',F)
        # print('d = ',d)

        ksi = np.linspace(-1,1,N,endpoint=False)

        for e in range(1,elements+1):
            x1 = he*(e-1)
            x2 = he*e
            xstep = (x2-x1)/N
            xtemp = np.zeros([N,1])
            for k in range(len(ksi)):
                xtemp[k] = x1 + k*xstep
            x = np.append(x,xtemp,axis=0)
            uhtemp = np.zeros([N,1])

            # set d1 and d2 for the given element
            d1 = d[e-1]
            if d1 == d[-1]:
                d2 = 0
            else:
                d2 = d[e]
            for k in range(len(ksi)):
                uhtemp[k] = 0.5*d1*(1-ksi[k]) + 0.5*d2*(1+ksi[k])

            uh = np.append(uh,uhtemp,axis=0)

        uh = np.delete(uh,0,0)
        x = np.delete(x,0,0)

        u = np.zeros([len(x),1])
        if i == 0:
            for j in range(len(x)):
                u[j] = 0.5*c**2 - 0.5*c*x[j]**2
        if i == 1:
            for j in range(len(x)):
                u[j] = (1/6.)*(1 - x[j]**3)
        if i == 2:
            for j in range(len(x)):
                u[j] = (1/12.)*(1 - x[j]**4)

        # Plot u(x) and uh(x) for the given combination of elements and load
        title = ('FEA solution for load case ' + fcase[i] + ' with '
            + str(elements) + ' elements')
        plt.figure()
        plt.title(title)
        plt.xlabel('Linear position along beam (x)')
        plt.ylabel('Displacement (u)')
        plt.plot(x,u,label='u(x)',linewidth=1,color='r')
        plt.plot(x,uh,label='uh(x)',linewidth=1,color='b',linestyle='--')
        plt.legend()
        plt.show()


        #######################################################################
        ## Part 2: Compute global error ##
        # Given constants for computing error integral with 3-point Gauss quad.
        ksii = [-np.sqrt(3./5.), 0, np.sqrt(3./5.)]
        wi = [5./9., 8./9., 5./9.]
        globerr = 0

        for e in range(1,elements+1):
            x1 = he*(e-1)
            x2 = he*e
            d1 = d[e-1]
            if d1 == d[-1]:
                d2 = 0
            else:
                d2 = d[e]
            xei = np.zeros([3,1])
            ue = np.zeros([3,1])
            uhe = np.zeros([3,1])
            for j in range(0,3):
                # Determine interpolated x values
                xei[j] = x1*0.5*(1 - ksii[j]) + x2*0.5*(1 + ksii[j])

            # Calculate the u value based on the appropriate x values
            if i == 0:
                for j in range(0,3):
                    ue[j] = 0.5*c**2 - 0.5*c*xei[j]**2
            if i == 1:
                for j in range(0,3):
                    ue[j] = (1/6.)*(1 - xei[j]**3)
            if i == 2:
                for j in range(0,3):
                    ue[j] = (1/12.)*(1 - xei[j]**4)

            # Calculate the uh values based on the appropriate ksi values
            for j in range(0,3):
                uhe[j] = 0.5*d1*(1-ksii[j]) + 0.5*d2*(1+ksii[j])

            # Perform Gauss quadrature for the current element:
            diff = np.zeros([3,1])
            for j in range(0,3):
                diff[j] = ue[j] - uhe[j]

            ab2 = np.zeros([3,1])
            for j in range(0,3):
                ab2[j] = (np.abs(diff[j]))**2

            gauss = np.zeros([3,1])
            for j in range(0,3):
                gauss[j] = ab2[j]*0.5*he*wi[j]

            # Sum up and add to global error
            total = 0
            for j in range(0,3):
                total += gauss[j]	# element total error
            globerr += total
        globerr = np.sqrt(globerr)
        print('global error = ',globerr)

        #######################################################################
        ## Part 3: Log-log plot ##
        if count == 0:
            evector[0] = globerr
        if count == 1:
            evector[1] = globerr
        if count == 2:
            evector[2] = globerr
        if count == 3:
            evector[3] = globerr
        count += 1

    convergence = 0
    convergence = (np.log10(evector[-1]/evector[0]))/(np.log10(n[-1]/n[0]))
    print("convergence = ",convergence)
    plt.figure()
    plt.title('Rate of convergence for load case ' + fcase[i])
    plt.xlabel('Number of elements (n)')
    plt.ylabel('Global error (e)')
    plt.text(200,0.0001,'Rate of convergence = '+str(convergence),horizontalalignment='center')
    plt.loglog(n,evector,linestyle='--',marker='o')
    plt.show()
