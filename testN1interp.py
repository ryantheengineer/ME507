# from Code1 import N1
import numpy as np

def N1(x1,x2,x):
    u = (x2 - x)/(x2 - x1)
    return u

x1 = 0
x2 = 1

# x = np.linspace(0,1,10,endpoint=True)
# u = np.zeros(len(x))

x = np.linspace(0,1,100)
print("x = ", x)
u = np.zeros(len(x))
print("u = ", u)

for i in range(len(x)):
    u[i] = N1(x1,x2,x[i])
print("u = ", u)
