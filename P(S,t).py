# %%
import pandas as pd 
import numpy as np 

left = 50
right = 100
dx = 1
#N = (right - left)/dx

T = 1
dt = 2 * 0.01 * 0.01
M = int(T/dt)

x = np.arange(left, right, dx)
x = np.append(x, right)
N = len(x)

r = 0.05 
sigma = 0.4
a = r - sigma * sigma/2
b = sigma * sigma/2
c = -1*r

w = np.zeros((M ) * (N )).reshape(M , N )
#w.shape
#w(0,:)
w[M-1,:] = (x - 50)*(100 - x)
for m in range(M-1,1,-1):
    print(m)
    t = m * dt
    w[m, 0] = 0 # Left  Boundary
    w[m, N-1] = 0 # Right Boundary
    for n in range(1,N-1):
        delta1 = -r*n/2 * (dt/dx) * (w[m,n+1]-w[m,n-1])
        delta2 = -b * n * n * (dt/dx/dx) * (w[m,n+1]-2*w[m,n]+w[m,n-1])
        delta3 = r * dt * w[m,n]
        delta = delta1 + delta2 + delta3 
        w[m-1,n] = w[m,n] + delta


# %%
import matplotlib.pyplot as plt 
P_S_0 = w[0,:]
P_S_1 = w[M-1,:]
plt.plot(x, P_S_0)
plt.plot(x, P_S_1)  
plt.legend(['V_x_0','V_x_1'])

# %%        
        
for i in range(10,0,-1):
    print(i)
        
        
        
        
        
        
        
        
        
        
        