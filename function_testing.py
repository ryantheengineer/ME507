import numpy as np
import matplotlib.pyplot as plt
from Code2 import Bap
# from Code2 import gaussquad
# from Code2 import fx
from Code2 import Bspline
from Code2 import knot
from Code2 import xAG
from Code2 import IEN
from Code2 import ID
from Code2 import LM

nel = 4
p = 3

s = knot(p,nel)
n = range(len(s)-p-1)
# print(A)

# Compute
xG = xAG(p,s,nel)
print('xG = ', xG)
xGlength = len(xG)
nodeindex = np.zeros([p+1,nel])

# for e in range(1,nel+1):
#     for a in range(1,p+2):
#         nodeindex[a-1,e-1] = LM(a,e,xGlength)
#
# print('nodeindex = ',nodeindex)
# The following section of code verifies the B-spline shape functions
plt.figure()
for e in range(1,nel+1):
    ksi=np.linspace(-1,1,10,endpoint=True)
    N = Bspline(e,p,nel,ksi)
    ksi+=2*(e-1)
    # print(ksi)
    for i in range(p+1):
        plt.plot(ksi,N[i,:])
    # plt.plot(ksi,N[0,:],ksi,N[1,:],ksi,N[2,:])
plt.show()



# p = [1,2,3]
# ksi = np.linspace(-1,1,100,endpoint=True)
# print('ksi = ',ksi)
# B = np.zeros([len(ksi),1])
# for order in p:
#     print('p = ',p)
#     plt.figure()
#     a = range(1,order+2)
#     print('a = ',a)
#     for curve in a:
#         for i in range(len(ksi)):
#             B[i] = Bap(curve,order,ksi[i])
#         plt.plot(ksi,B)
#     plt.show()
# ## Bap functions as it's supposed to

# def fx(x):
#     f = np.zeros([len(x),1])
#     # print(f)
#     for i in range(len(x)):
#         # print('x[i] = ',x[i])
#         # print('x[i]**2 = ',x[i]**2)
#         f[i] = x[i]**2
#     return f
#
# x = np.linspace(-1,1,100,endpoint=True)
# g = fx(x)
# plt.figure()
# plt.plot(x,g)
# plt.show()
# # print(g)
#
# nint = 3
# areabelowcurve = gaussquad(nint)
# print('areabelowcurve = ',areabelowcurve)
