import numpy as np
import matplotlib.pyplot as plt
# from Code2 import Bap
# from Code2 import gaussquad
# from Code2 import fx
from Code2 import knot

nel = 10
p = 1
knot(p,nel)









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
