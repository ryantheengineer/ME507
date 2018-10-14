# from Code1 import N1
import numpy as np
import matplotlib.pyplot as plt

def N1(x1,x2,x):
    u = (x2 - x)/(x2 - x1)
    return u

def N2(x1,x2,x):
    u = (x - x1)/(x2 - x1)
    return u

n = 4
x1 = 0
x2 = 1
d = [0.3,0.15,0.1,0.075]

# x = np.linspace(0,1,10,endpoint=True)
# u = np.zeros(len(x))

x = np.linspace(0,1,100)
# print("x = ", x)
u = np.zeros(len(x))
# print("u = ", u)

# for i in range(len(x)):
#     u[i] = N1(x1,x2,x[i])
# print("u = ", u)

he = 1/float(n)
print("he = %f") % he
uh = np.zeros(len(x))

for el in range(1,n+1):
    x1 = (el-1)*he
    x2 = el*he
    print("x1 = %f") % x1
    print("x2 = %f") % x2
    d1 = d[el-1]
    if d1 == d[-1]:
        d2 = 0
    else:
        d2 = d[el]
    print(type(d1))
    print("d1 = ", d1)
    print("d2 = ", d2)


    for j in range(len(x)):
        # print(type(j))
        if x[j] >= x1 and x[j] < x2:
            print(x[j])
            uh[j] += d1*N1(x1,x2,x[j]) + d2*N2(x1,x2,x[j])
        else:
            uh[j] += 0

print(uh)
plt.figure()
plt.plot(x,uh)
plt.show()
