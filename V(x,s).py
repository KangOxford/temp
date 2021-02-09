# %%
import pandas as pd 
import numpy as np 

left = np.log(50)
right = np.log(100)
dx = 1/100
#N = (right - left)/dx

T = 1
k = 1/4
dt = k * dx * dx
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
w[0,:] = (np.exp(x) - 50)*(100 - np.exp(x))
for m in range(M-1):
    t = m * dt
    w[m, 0] = 0 # Left  Boundary
    w[m, N-1] = 0 # Right Boundary
    for n in range(1,N-1):
        delta1 = a/2 * (dt/dx) * (w[m,n+1]-w[m,n-1])
        delta2 = b * (dt/dx/dx) * (w[m,n+1]-2*w[m,n]+w[m,n-1])
        delta3 = c * dt * w[m,n]
        delta = delta1 + delta2 + delta3 
        w[m+1,n] = w[m,n] + delta

print(w[0,1:])
# %%
# import matplotlib.pyplot as plt 
# V_x_0 = w[0,:]
# V_x_1 = w[M-1,:]
# plt.plot(x, V_x_0)
# plt.plot(x, V_x_1)  
# plt.legend(['V_x_0','V_x_1'])

# %%        
        
        
        
        
        
        
        
        
        
        
        
        
        