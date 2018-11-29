'''Code2.py: FEA code for ME 507'''

import numpy as np
import matplotlib.pyplot as plt
import sys
import math

## FUNCTIONS: ##
# Bezier curves
def Bap(a,p,ksi):
    num = math.factorial(p)
    den = math.factorial(a-1)*math.factorial(p+1-a)
    pa_1 = num/den
    B = (1./(2.**p))*pa_1*((1-ksi)**(p-(a-1)))*((1+ksi)**(a-1))
    return B

# 1st derivative of Bap with respect to ksi (as given by Wolfram Alpha)
def Bap1(a,p,ksi):
    num1 = (a-1.)*math.factorial(p)*((ksi+1.)**(a-2.))*((1.-ksi)**(-a+p+1.))
    den1 = (2.**p)*math.factorial(a-1)*math.factorial(-a+p+1.)
    num2 = math.factorial(p)*(-a+p+1.)*((ksi+1.)**(a-1.))*((1.-ksi)**(p-a))
    den2 = den1
    B1 = (num1/den1) - (num2/den2)
    return B1

# 2nd derivative of Bap with respect to ksi (as given by Wolfram Alpha)
def Bap2(a,p,ksi):
    term1 = 1./(math.factorial(a - 1.)*math.factorial(-a + p + 1.))
    term2 = (a-2.)*(a-1.)*(2.**(-p))*((ksi+1.)**(a-3.))*((1.-ksi)**(-a+p+1.))
    term3 = (a-1.)*(2.**(1.-p))*(-a+p+1.)*((ksi+1.)**(a-2.))*((1.-ksi)**(p-a))
    term4 = (2.**(-p))*(p-a)*(-a+p+1.)*((ksi+1.)**(a-1.))*((1.-ksi)**(-a+p-1.))
    B2 = term1*math.factorial(p)*(term2 - term3 + term4)
    return B2

# Extraction operators, which we use to convert Bezier curves to B-splines
# Define Ce for p = 2
def Ce2(e,nel): # NOTE: may need to check indexing since elements need to be adjusted to start at 0?
    if e == 1:
        Ce = np.array([[1., 0, 0],
                       [0, 1., 0.5],
                       [0, 0, 0.5]])
    if e > 1 and e < nel:
        Ce = np.array([[0.5, 0, 0],
                       [0.5, 1., 0.5],
                       [0, 0, 0.5]])
    if e == nel:
        Ce = np.array([[0.5, 0, 0],
                       [0.5, 1., 0],
                       [0, 0, 1.]])
    return Ce

# Define Ce for p = 3
def Ce3(e,nel):
    if e == 1:
        Ce = np.array([[1., 0, 0, 0],
                       [0, 1., 0.5, 0.25],
                       [0, 0, 0.5, 7./12.],
                       [0, 0, 0, 1./6.]])
    if e == 2:
        Ce = np.array([[0.25, 0, 0, 0],
                       [7./12., 2./3., 1./3., 1./6.],
                       [1./6., 1./3., 2./3., 2./3.],
                       [0, 0, 0, 1./6.]])
    if e > 2 and e < nel-1:
       Ce = np.array([[1./6., 0, 0, 0],
                      [2./3., 2./3., 1./3., 1./6.],
                      [1./6., 1./3., 2./3., 2./3.],
                      [0, 0, 0, 1./6.]])
    if e == nel-1:
        Ce = np.array([[1./6., 0, 0, 0],
                       [2./3., 2./3., 1./3., 1./6.],
                       [1./6., 1./3., 2./3., 7./12.],
                       [0, 0, 0, 0.25]])
    if e == nel:
        Ce = np.array([[1./6., 0, 0, 0],
                       [7./12., 0.5, 0, 0],
                       [0.25, 0.5, 1., 0],
                       [0, 0, 0, 1.]])
    return Ce

def fx(x):  # x must be a vector
    f = np.zeros([len(x),1])
    for i in range(len(x)):
        f[i] = x[i]**2
    return f

# Set up quadrature rule (not integration): ksi, w
def gaussksi(nint):
    ksiint = np.zeros([nint,1])
    if nint == 1:
        ksiint[0] = 0
    elif nint == 2:
        ksiint[0] = -1./np.sqrt(3.)
        ksiint[1] = 1./np.sqrt(3.)
    elif nint == 3:
        ksiint[0] = -np.sqrt(3./5.)
        ksiint[1] = 0
        ksiint[2] = np.sqrt(3./5.)
    return ksiint

def gaussW(nint):
    W = np.zeros([nint,1])
    if nint == 1:
        W[0] = 0
    elif nint == 2:
        W[0] = 1
        W[1] = 1
    elif nint == 3:
        W[0] = 5./9.
        W[1] = 8./9.
        W[2] = 5./9.
    return W


# Set up LM, ID, IEN arrays
# a must be a column vector containing 1,2,..., # of shape functions
def IEN(a,e):
    A = 0
    if a == 1:
        A = e
    if a == 2:
        A = e + 1
    if a == 3:
        A = e + 2
    if a == 4:
        A = e + 3
    # if a == 5:
    #     A = e + 4
    return A

# Take in global node A. If the global node corresponds to an inactive node, in
# the case of this problem, at L = 1, then the ID array outputs 0 and continues
# numbering afterward. So maybe have it take in a list of inactive node numbers?
def ID(node,xGlength):
    eq = 0
    if node == xGlength or node == xGlength-1:  # NOTE: This is a change from the original Code2. Check on this.
        eq = 0
    else:
        eq = node
    return eq

def LM(a,e,xGlength):
    A = IEN(a,e)
    P = ID(A,xGlength)
    return P

# Compute node locations
# Create knot vector
def knot(p,nel):
    he = 1./nel
    s = np.zeros([2*p+nel+1,1])
    temp = 0
    for i in range(len(s)):
        if i <= p or i > nel+p:
            temp += 0
        else:
            temp += he
        s[i] = temp
    return s

# Use knot vector and p value to determine node locations
def xAG(p,s,nel):
    n = range(len(s)-p-1)
    xG = np.zeros([len(n),1])
    for A in n:
        # print('A = ',A)
        for i in range(A+1,A+p+1):
            # print('s[i] = ',s[i])
            xG[A] += (1.0/p)*s[i]
    return xG


# Define Bsplines for a given element
def Bspline(e,p,nel,ksi):
    Be = np.zeros([p+1,len(ksi)])
    for a in range(1,p+2):
        for i in range(len(ksi)):
            Be[a-1,i] = Bap(a,p,ksi[i])
            # This gives us our len(ksi)x3 or len(ksi)x4 arrays of Bezier curve
            # interpolation values in the parent domain. Next we need to convert
            # these into B-splines
    if p == 1:
        Ce = np.array([[1, 0],
                       [0, 1]])
    if p == 2:
        Ce = Ce2(e,nel)
    if p == 3:
        Ce = Ce3(e,nel)

    # Ne should be an array the same size as Be
    Ne = np.zeros([p+1,len(ksi)])
    # Create vectors of B-spline values, depending on the ksi value on the parent domain
    Becol = np.zeros([p+1,1])
    Necol = np.zeros([p+1,1])
    for i in range(len(ksi)):
        for j in range(p+1):
            Becol[j] = Be[j,i]
        # print(np.shape(Becol))
        # print(np.shape(Ce*Becol))
        Necol = np.matmul(Ce,Becol)
        # print(np.shape(Necol))
        for j in range(p+1):
            Ne[j,i] = Necol[j]
    return Ne

################################ BEGIN MAIN CODE ###############################

if __name__ == "__main__":
    ######### INPUT ###########
    # nel = [1, 10, 100]
    nel = [1]
    # p = [2, 3]
    p = [3]
    E = 1000000.0
    b = 0.005
    h = 0.005
    I = (b*h**3)/12.
    f = 10*h**3

    ######### SETUP ###########
    # Set up gauss quadrature rule
    nint = 3
    ksiint = gaussksi(nint)
    W = gaussW(nint)

    for P in p:
        # print('\n\n')
        # print('p = ',P)
        evector = np.zeros([len(nel),1])
        hvector = np.zeros([len(nel),1])
        for i in range(len(hvector)):
            hvector[i] = 1./nel[i]
        nodevector = np.zeros([len(nel),1])
        count = 0

        for elements in nel:
            he = 1./elements
            nen = P + 1     # number of element nodes
            knotvector = knot(P,elements)
            Narray = np.zeros([P+1,elements*nint])
            N1array = np.zeros([P+1,elements*nint])
            N2array = np.zeros([P+1,elements*nint])
            # dNdx = np.zeros([P+1,elements*nint])
            # d2Ndx2 = np.zeros([P+1,elements*nint])
            # print('Narray dimensions: ',np.shape(Narray))
            xarray = np.zeros([1])
            x1array = np.zeros([1])
            x2array = np.zeros([1])
            # print('knotvector = ',knotvector)
            # xG = []
            xG = xAG(P,knotvector,elements)  # nodes
            # print('xG = ',xG)
            xGlength = len(xG)
            nodevector[count] = xGlength
            activenodes = len(xG)-1
            # print('xG length: ',len(xG))
            print('# of active nodes: ',activenodes)
            K = np.zeros([activenodes,activenodes])
            F = np.zeros([activenodes,1])
            col = 0
            for e in range(1,elements+1):
                print('\n')
                print('Element: ',e-1)
                if P == 2:
                    Ce = Ce2(e,elements)
                elif P == 3:
                    Ce = Ce3(e,elements)
                # print('Ce = ',Ce)
                if elements == 1 and P == 2:
                    Ce = np.array([[1., 0., 0.],
                                   [0., 1., 0.],
                                   [0., 0., 1.]])
                if elements < 5 and P == 3:
                    Ce = np.array([[1.,0.,0.,0.],
                                   [0.,1.,0.,0.],
                                   [0.,0.,1.,0.],
                                   [0.,0.,0.,1.]])

                fe = np.zeros([nen,1])
                ke = np.zeros([nen,nen])

                # INTEGRATION LOOP
                for i in range(1,nint+1):
                    print('\n')
                    print('\ti = ' + str(i-1))
                    wi = W[i-1]
                    Be = np.zeros([P+1,1])
                    B1e = np.zeros([P+1,1])
                    B2e = np.zeros([P+1,1])
                    x = 0.0
                    x1 = 0.0
                    x2 = 0.0
                    # iterate on the Bernstein polynomials for the element
                    for a in range(1,P+2):
                        Be[a-1] = Bap(a,P,ksiint[i-1])
                        B1e[a-1] = Bap1(a,P,ksiint[i-1])
                        B2e[a-1] = Bap2(a,P,ksiint[i-1])
                    Ne = np.matmul(Ce,Be)
                    N1e = np.matmul(Ce,B1e) # 1st derivative of Ne in parent domain
                    N2e = np.matmul(Ce,B2e) # 2nd derivative of Ne in parent domain
                    # print('xG = ', xG[a-1])
                    print('Ne = ' + str(Ne))
                    print('N1e = ' + str(N1e))

                    # insert Ne into Narray
                    for row in range(P+1):
                        # column of Narray is the integration point
                        # row of Narray is the shape function number
                        Narray[row,col] = Ne[row] # Narray is in the parent domain, and includes all element values put into one array for the whole beam
                        # print('Narray = ',Narray)
                    for row in range(P+1):
                        N1array[row,col] = N1e[row]
                    for row in range(P+1):
                        N2array[row,col] = N2e[row]


                    for a in range(1,P+2):
                        x += xG[a-1+(e-1)]*Ne[a-1]
                    xarray = np.append(xarray,x) # x based on ksi; includes all values, not just a single element
                    print('x = ' + str(x))
                    # print('xarray = ',xarray)

                    for a in range(1,P+2):
                        x1 += xG[a-1+(e-1)]*N1e[a-1]
                    x1array = np.append(x1array,x1)

                    for a in range(1,P+2):
                        x2 += xG[a-1+(e-1)]*N2e[a-1]
                    x2array = np.append(x2array,x2)

                    # QUESTION: for loop to create dN/dx? Is dN/dx an array or a specific value?
                    dNdx = N1e*(x1)**(-1)

                    d2Ndx2 = (N2e - N1e*x2)*((x1**2)**(-1))

                    # fi = f
                    fi = x**2
                    print('f(x) = ' + str(fi))

                    # calculate fe vector
                    for a in range(0,nen):
                        fe[a] += Ne[a]*fi*(he/2.)*wi # QUESTION: is this defined correctly?
                        for b in range(0,nen):
                            ke[a,b] += d2Ndx2[a]*E*I*d2Ndx2[b]*((2./he)**3)*wi # FIXME: This is almost correct, but not quite

                    col += 1
                    print('fe = ' + str(fe))
                    print('ke = ' + str(ke))
                    # print('\n')

                for a in range(1,nen+1):
                    # print('a = ',a)
                    # print('LM(a,e)= ',LM(a,e,xGlength))
                    if LM(a,e,xGlength) > 0:
                        F[LM(a,e,xGlength)-1] += fe[a-1]
                        # print('F = ',F)
                        for b in range(1,nen+1):
                            if LM(b,e,xGlength) > 0:
                                K[LM(a,e,xGlength)-1,LM(b,e,xGlength)-1] += ke[a-1,b-1]

            xarray = np.delete(xarray,0)
            x1array = np.delete(x1array,0)
            x2array = np.delete(x2array,0)
            # print('xarray = ',xarray)
            # print('x1array = ',x1array)
            # print('x2array = ',x2array)
            # print('Narray = ',Narray)
            print('F = ',F)
            print('K = ',K)

            d = np.zeros([activenodes,1])
            d = np.linalg.solve(K,F)

            # append 0 on the end for calculating uh
            d = np.append(d,[0.,0.])
            # d = np.append(d,0.)
            print('d = ',d)

            # Calculate uh(x)
            uh = np.zeros([len(xarray),1])
            for x in range(len(uh)):
                loc = x/nint # integer devision returns the element number
                for A in range(P+1):
                    uh[x] += d[loc+A]*Narray[A,x]

            # Calculate u(x)
            u = np.zeros([len(xarray),1])
            for j in range(len(xarray)):
                u[j] = ((f*xarray[j]**2)/(24.*E*I))*(2.+(2.-xarray[j])**2)

            # Plot u(x) and uh(x) for the given combination of elements and load
            title = ('FEA solution for p = ' + str(P) + ' with '
                + str(elements) + ' elements')
            plt.figure()
            plt.title(title)
            plt.xlabel('Linear position along beam (x)')
            plt.ylabel('Displacement (u)')
            plt.plot(xarray,u,label='u(x)',linewidth=1,color='r')
            plt.plot(xarray,uh,label='uh(x)',linewidth=1,color='b',linestyle='--')
            plt.legend()
            plt.show()

            # Calculate the global error
            i = 0
            total = 0.0
            globerr = 0.0
            diff = np.zeros([len(uh),1])
            ab2 = np.zeros([len(uh),1])
            gauss = np.zeros([len(uh),1])
            for x in range(len(uh)):
                diff[x] = u[x] - uh[x]
                ab2[x] = (np.abs(diff[x]))**2
                wi = W[i]   # want indices 0 through 2, repeating 5 times
                gauss[x] = ab2[x]*0.5*he*wi
                globerr += gauss[x]
                i += 1  # increase the integration point index on the element level
                if i == 3:
                    i = 0
            globerr = np.sqrt(globerr)
            evector[count] = globerr
            count += 1

        convergence = 0
        convergence = (np.log10(evector[-1]/evector[0]))/(np.log10(hvector[-1]/hvector[0]))

        ## Log-log plots ##
        plt.figure()
        plt.title('Rate of convergence for p = ' + str(P))
        plt.xlabel('Element size (h)')
        plt.ylabel('Global error (e)')
        # plt.text(hvector[1],evector[2],'Rate of convergence = '+str(convergence),horizontalalignment='center')
        plt.loglog(hvector,evector,linestyle='--',marker='o')
        plt.show()

        plt.figure()
        plt.title('Error vs. Nodes for p = ' + str(P))
        plt.xlabel('Number of nodes')
        plt.ylabel('Global error (e)')
        plt.loglog(nodevector,evector,linestyle='--',marker='o')
        plt.show()
