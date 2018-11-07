'''Code2.py: FEA code for ME 507'''

import numpy as np
import matplotlib.pyplot as plt
import sys
import math

## INPUT: ##
# P (What does this mean exactly? Is it referring to the P that comes out of the LM function?)

## SETUP: ##
# Bezier curves
def Bap(a,p,ksi):
    num = math.factorial(p)
    den = math.factorial(a-1)*math.factorial(p+1-a)
    pa_1 = num/den
    B = (1./2.**p)*pa_1*((1-ksi)**(p-(a-1)))*((1+ksi)**(a-1))
    return B

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

# def fx(x):  # x must be a vector
#     f = np.zeros([len(x),1])
#     for i in range(len(x)):
#         f[i] = x[i]**2
#     return f

# Set up quadrature rule: ksi, w
# NOTE: Gauss quadrature might need to happen on the global x coordinate rather
# than the parent domain. Come back to this later
# def gaussquad(nint):
#     W = np.zeros([nint,1])
#     ksiint = np.zeros([nint,1])
#     if nint == 1:
#         ksiint[0] = 0
#         W[0] = 0
#     elif nint == 2:
#         ksiint[0] = -1./np.sqrt(3.)
#         ksiint[1] = 1./np.sqrt(3.)
#         W[0] = 1
#         W[1] = 1
#     elif nint == 3:
#         ksiint[0] = -np.sqrt(3./5.)
#         ksiint[1] = 0
#         ksiint[2] = np.sqrt(3./5.)
#         W[0] = 5./9.
#         W[1] = 8./9.
#         W[2] = 5./9.
#
#     integration = 0
#     for i in range(nint):
#         print('integration = ',integration)
#         integration += fx(ksiint[i])*W[i]
#         print('integration = ',integration)
#     return integration  # This gives the integral only in the parent domain ([-1,1]) and not in the global domain


# Set up LM, ID, IEN arrays
# a must be a 2x1 vector containing 1,2
def IEN(a,e):
    if a == 1:
        A = e
    if a == 2:
        A = e + 1
    return A

def ID(a,e,n):
    if e == n and a == 2:
        A = 0
    return A

def LM(a,e,n):
    P = ID(IEN(a,e),a,e,n)
    return P

# Compute node locations
# Create knot vector
def knot(p,nel):
    he = 1./nel
    # print('he = ',he)
    # he = 1.0
    s = np.zeros([p+nel+1+p,1]) # knot vector should have length of 2*p + # of nodes
    temp = 0
    for i in range(len(s)):
        if i <= p or i > nel+p:
            temp += 0
        else:
            temp += he
        s[i] = temp

    # print('s = ', s)
    return s

def xAG(p,s):
    n = range(len(s)-p-1)
    xG = np.zeros([len(n),1])
    for A in n:
        for i in range(A+1,A+p+1):
            xG[A] += (1.0/p)*s[i]
    return xG


# Define Bsplines for a given element NOTE: MAY NEED TO CHANGE SO YOU FEED IN A SINGLE VALUE AND GET THE SINGLE N INTERPOLATED VALUE OUT FOR EACH N
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


nel = [1,10,100,1000]



## NEED TO DEFINE B-SPLINE INTERPOLATION FUNCTIONS HERE

# nint should be p + 1. This says how many integration points per element
nint = 3    # NOTE: check this value and make sure it isn't supposed to change







################################ PSEUDO-CODE: ##################################
# Initialize K and F
# for elements in nel:
#     K = np.zeros([elements,elements])
#     F = np.zeros([elements,1])
# # Define f(x)
#     # fx = x**2
# # for e = 1,2,...,nel:
#     for e in range(0,elements):
#         # initialize ke and fe
#         # for i = 1,2,...,nint (where nint is the # of integration points): (THIS IS THE INTEGRATION LOOP)
#         for i in range(0,nint):
#             # evaluate B[a][p](ksi[i])
#             # B[a][p] = Bap(a,p,ksi)  # FIXME: Not sure if this is right
#             # compute N = C[e]*B(ksi[i])
#             # evaluate derivative of B(ksi[i]) with respect to ksi
#             # Take derivative of N with respect to ksi, which is equal to C[e]*derivative of B with respect to ksi
#             # compute
#
#
#
#
#
#         # Assemble the element matrices into the global matrices
#         for a in range(0,nel):  # NOTE: In the notes there are two different values, nel and nen (are these the same thing?)
#             if LM(a,e,elements) > 0:
#                 F[LM(a,e,elements)] = F[LM(a,e,elements)] + fe[a]
#             for b in range(0,nel):
#                 if LM(b,e,elements) > 0:
#                     K[LM(a,e,elements),LM(b,e,elements)] += ke[a,b] # NOTE: not sure on the syntax, but the point is to add the ke entry into K at the right position
#
#
#     # end of e loop
#
#     # Solve Kd = F
