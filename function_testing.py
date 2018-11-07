import numpy as np
import matplotlib.pyplot as plt
from Code2 import Bap
# from Code2 import gaussquad
# from Code2 import fx
from Code2 import Bspline
from Code2 import knot

nel = 10
p = 2

s = knot(p,nel)
n = range(len(s)-p-1)
# print(A)

# Graville Abscissae
xG = np.zeros([len(n),1])
for A in n:
    for i in range(A+1,A+p+1):
        xG[A] += (1.0/p)*s[i]
print('xG = ',xG) # It looks like this snippet works properly now

# xs = np.zeros([len(xG),1])
# for knot in s:
#     for A in n:
#         xs[knot] += xG[A]*


# if n = len(s), then there are n-p-1 shape functions
# Convert to B-spline:
# for each element:
#     assemble vector of

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
