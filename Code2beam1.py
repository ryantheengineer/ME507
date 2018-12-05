'''Code2beam1.py: FEA code for ME 507
   Part 1 of the beam bending problem'''

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
def Bap2(a,p,x):
    B2 = (((1. - x)**(-1. - a + p))*((1. + x)**(-3. + a))*(4. + 4.*(a**2.)
        - 4.*x + (p**2.)*(1. + x)**2. - 4.*a*(2. + p - x + p*x)
        + p*(3. + 2.*x - x**2.))*math.factorial(p))/((2.**p)*math.factorial(-1. + a)*math.factorial(1. - a + p))

    return B2


# Extraction operators, which we use to convert Bezier curves to B-splines
# Define Ce for p = 2
def Ce2(e,nel):
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
    if node == xGlength or node == xGlength-1:
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
        Necol = np.matmul(Ce,Becol)
        for j in range(p+1):
            Ne[j,i] = Necol[j]
    return Ne

################################ BEGIN MAIN CODE ###############################

if __name__ == "__main__":
    ######### INPUT ###########
    nel = [1, 10, 100]
    p = [2, 3]
    E = 1000000.0
    b = 0.005
    h = 0.005
    f = 10.*(h**3.)
    I = (b*(h**3.))/12.

    ######### SETUP ###########
    # Set up gauss quadrature rule
    nint = 3
    ksiint = gaussksi(nint)
    W = gaussW(nint)

    # Set up holding array for tip deflections. Rows for p and columns for nel
    deflections = np.zeros([len(p),len(nel)])
    exactsol = np.zeros([len(nel),1])

    for elements in nel:
        # plt.figure()
        print('\nn = ' + str(elements))

        for P in p:
            print('\n\tp = ' + str(P))
            evector = np.zeros([len(nel),1])
            hvector = np.zeros([len(nel),1])
            for i in range(len(hvector)):
                hvector[i] = 1./nel[i]
            nodevector = np.zeros([len(nel),1])
            count = 0

            he = 1./elements
            nen = P + 1     # number of element nodes
            knotvector = knot(P,elements)
            Narray = np.zeros([P+1,elements*nint])
            N1array = np.zeros([P+1,elements*nint])
            N2array = np.zeros([P+1,elements*nint])

            xarray = np.zeros([1])
            x1array = np.zeros([1])
            x2array = np.zeros([1])

            xG = xAG(P,knotvector,elements)  # nodes
            xGlength = len(xG)
            nodevector[count] = xGlength
            activenodes = len(xG)-1

            K = np.zeros([activenodes-1,activenodes-1])
            F = np.zeros([activenodes-1,1])
            col = 0

            for e in range(1,elements+1):
                if P == 2:
                    Ce = Ce2(e,elements)
                if P == 3:
                    Ce = Ce3(e,elements)

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

                    # insert Ne into Narray
                    for row in range(P+1):
                        # column of Narray is the integration point
                        # row of Narray is the shape function number
                        Narray[row,col] = Ne[row]
                    for row in range(P+1):
                        N1array[row,col] = N1e[row]
                    for row in range(P+1):
                        N2array[row,col] = N2e[row]
                    col += 1

                    for a in range(1,P+2):
                        x += xG[a-1+(e-1)]*Ne[a-1]
                    xarray = np.append(xarray,x) # x based on ksi; includes all values, not just a single element

                    for a in range(1,P+2):
                        x1 += xG[a-1+(e-1)]*N1e[a-1]
                    x1array = np.append(x1array,x1)

                    for a in range(1,P+2):
                        x2 += xG[a-1+(e-1)]*N2e[a-1]
                    x2array = np.append(x2array,x2)

                    fi = f

                    # calculate fe vector
                    for a in range(0,nen):
                        fe[a] += Ne[a]*fi*(he/2.)*wi
                        for b in range(0,nen):
                            ke[a,b] += N2e[a]*E*I*N2e[b]*((2./he)**3.)*wi

                for a in range(1,nen+1):
                    if LM(a,e,xGlength) > 0:
                        F[LM(a,e,xGlength)-1] += fe[a-1]
                        for b in range(1,nen+1):
                            if LM(b,e,xGlength) > 0:
                                K[LM(a,e,xGlength)-1,LM(b,e,xGlength)-1] += ke[a-1,b-1]

            xarray = np.delete(xarray,0)
            xushift = np.zeros([len(xarray),1])
            for x in range(len(xushift)):
                xushift[x] = -xarray[x] + 1.
            x1array = np.delete(x1array,0)
            x2array = np.delete(x2array,0)
            # print('xarray = ',xarray)
            # print('x1array = ',x1array)
            # print('x2array = ',x2array)
            # print('Narray = ',Narray)
            # print('K = ',K)
            # print('F = ',F)

            d = np.zeros([activenodes,1])
            d = np.linalg.solve(K,F)
            # append 0 on the end for calculating uh
            d = np.append(d,[0.,0.])
            # print('d = ',d)

            # Calculate uh(x)
            uh = np.zeros([len(xarray),1])
            for x in range(len(uh)):
                loc = x/nint # integer devision returns the element number
                for A in range(P+1):
                    uh[x] += d[loc+A]*Narray[A,x]

            maxtip = uh[0]
            print('\tMax tip deflection = ' + str(maxtip))

            deflections[p.index(P),nel.index(elements)] = maxtip
            count += 1

        exactsol[nel.index(elements)] = f/(8.0*E*I)

    # Plot the results asked for in part 2, problem 1:
    plt.title('Max tip deflection calculations')
    plt.xlabel('Number of elements (n)')
    plt.ylabel('Tip deflection (u)')
    plt.plot(nel,deflections[0,:],label='p = 2',linewidth=1,color='b',
        linestyle='--',marker='o')
    plt.plot(nel,deflections[1,:],label='p = 3',linewidth=1,color='g',
        linestyle='--',marker='o')
    plt.plot(nel,exactsol,label='Exact',linewidth=1,color='r')
    plt.legend()
    plt.show()
