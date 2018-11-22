'''Code2.py: FEA code for ME 507'''

import numpy as np
import matplotlib.pyplot as plt
import sys
import math
# import math.factorial as fac

## INPUT: ##
# P (What does this mean exactly? Is it referring to the P that comes out of the LM function?)

## SETUP: ##
# Bezier curves
def Bap(a,p,ksi):
    num = math.factorial(p)
    den = math.factorial(a-1)*math.factorial(p+1-a)
    pa_1 = num/den
    B = (1./(2.**p))*pa_1*((1-ksi)**(p-(a-1)))*((1+ksi)**(a-1))
    return B

# 1st derivative of Bap with respect to ksi (as given by Wolfram Alpha)
def Bap1st(a,p,ksi):
    num1 = (a-1)*math.factorial(p)*((ksi+1)**(a-2))*((1-ksi)**(-a+p+1))
    den1 = 2*p*math.factorial(a-1)*math.factorial(-a+p+1)
    num2 = math.factorial(p)*(-a+p+1)*((ksi+1)**(a-1))*((1-ksi)**(p-a))
    den2 = den1
    B1 = (num1/den1) - (num2/den2)
    return B1

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
    if a == 1:
        A = e
    if a == 2:
        A = e + 1
    if a == 3:
        A = e + 2
    return A

# Take in global node A. If the global node corresponds to an inactive node, in
# the case of this problem, at L = 1, then the ID array outputs 0 and continues
# numbering afterward. So maybe have it take in a list of inactive node numbers?
def ID(node,xGlength):
    if node == xGlength:
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
def knot(p,nel):    # FIXME: It looks like this might be creating an incorrect length, resulting in a non-1 end node
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
def Bspline(e,p,nel,ksi):   # give it a ksi vector like the one below
    # ksi = np.linspace(-1,1,100,endpoint=True):
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



if __name__ == "__main__":
    ###########################
    ######### INPUT ###########
    ###########################
    # nel = [1, 10, 100, 1000]  # FIXME: How do you deal with nel = 1?
    # nel = [10, 100, 1000]
    # nel = [10]
    nel = [3]
    # p = [2, 3]
    p = [2]
    a = np.array([[1],[2]])

    ###########################
    ######### SETUP ###########
    ###########################

    # Set up gauss quadrature rule
    nint = 3 # nint = p+1 is normally sufficient. We don't have the nint = 4 gauss rule so we will use nint = 3
    ksiint = gaussksi(nint)
    W = gaussW(nint)
    # print('ksiint = ', ksiint)
    # print('W = ', W)

    # NOTE: need to complete setup of Bap, Ce, quadrature rule, setup arrays, and node locations (including knot vectors)

    # Set up LM array
    # for elements in nel:
    # LMarray = LM()

    # Integration loop
    for elements in nel:

        # x = np.linspace(0,1,10*elements,endpoint=True)
        # f = fx(x)
        he = 1./elements

        for P in p:
            print('\n\n')
            print('p = ',P)
            nen = P + 1     # number of element nodes
            knotvector = knot(P,elements)
            # print('knotvector = ',knotvector)
            # xG = []
            xG = xAG(P,knotvector,elements)  # nodes
            print('xG = ',xG)
            xGlength = len(xG)
            activenodes = len(xG)-1
            print('xG length: ',len(xG))
            print('# of active nodes: ',activenodes)
            K = np.zeros([activenodes,activenodes])
            F = np.zeros([activenodes,1])

            for e in range(1,elements+1):
                print('\n')
                print('Element: ',e-1)
                if P == 2:
                    Ce = Ce2(e,elements)
                elif P == 3:
                    Ce = Ce3(e,elements)
                # print('Ce = ',Ce)
                if elements == 1:
                    Ce = np.array([[1, 0],[0, 1]])
                # ke = np.zeros([2,2])    # this size is incorrect
                # fe = np.zeros([2,1])    # this size is incorrect
                fe = np.zeros([nen,1])  # Check these dimensions later
                ke = np.zeros([nen,nen])

                # iterate on each integration point
                for i in range(1,nint+1):
                    wi = W[i-1]
                    Be = np.zeros([P+1,1])
                    B1e = np.zeros([P+1,1])
                    x = 0.0
                    # print('\n')
                    # print('Integration point (i): ',i-1)
                    # iterate on the Bernstein polynomials for the element
                    for a in range(1,P+2):
                        Be[a-1] = Bap(a,P,ksiint[i-1])
                        B1e[a-1] = Bap1st(a,P,ksiint[i-1])
                        Ne = np.matmul(Ce,Be)
                        Ne1 = np.matmul(Ce,B1e) # 1st derivative of Ne
                        # print('xG = ', xG[a-1])
                        # print('Ne = ', Ne[a-1])

                    # print('B array: ',Be)
                    # print('dB_dxi array: ',B1e)
                    # print('\n')
                    # print('N array: ',Ne)
                    # print('dN_dxi array: ',Ne1)
                    for a in range(1,P+2):
                        x += xG[a-1+(e-1)]*Ne[a-1]
                    # print('x = %f') % x
                    fi = x**2
                    # print('f(x) = ',fi)

                    # calculate fe vector

                    for a in range(0,nen):
                        fe[a] += Ne[a]*fi*(he/2.)*wi
                        for b in range(0,nen):
                            ke[a,b] += Ne1[a]*Ne1[b]*(2./he)*wi

                    # print('fe = ',fe)
                    # print('\n')
                    # print('ke = ',ke)

                for a in range(1,nen+1):
                    # print('LM(a,e)= ',LM(a,e,xGlength))
                    if LM(a,e,xGlength) > 0:
                        F[LM(a,e,xGlength)-1] += fe[a-1]
                        # print('F is currently: ',F)
                        for b in range(1,nen+1):
                            if LM(b,e,xGlength) > 0:
                                print('LM(a,e) = ',LM(a,e,xGlength))
                                print('LM(b,e) = ',LM(b,e,xGlength))
                                K[LM(a,e,xGlength)-1,LM(b,e,xGlength)-1] += ke[a-1,b-1]
                                # print('K is currently: ',K)

        d = np.linalg.solve(K,F)
        print('d = ',d)
