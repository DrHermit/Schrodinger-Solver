#Sum_i,j,k Tx[i’,i] x[i,j,k] + Ty[j’,j]x[i,j,k] + Tz[k’,k]x[i,j,k] = y[i’,j’,k’]
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
from scipy.sparse.linalg import LinearOperator
import math

#parameters of variables
#number of grid points
N = 201
#mass
m = 1
xmax = 5
xmin = -5
L = xmax - xmin
x = np.linspace(xmin, xmax, N)
#for N even
C = (0.5*(math.pi*N/L)**2)/m
#number of eigenvalues
n_eigs = 5

#define potential function
V = np.linspace(xmin,xmax,N)
for i in range(N):
    V[i] = .5*V[i]**2

H = sparse.eye(N, N, format='lil')*0

#Defining Kinetic Operator
if N%2 == 0:
    for i in range(N):
        for j in range(N):
	        if j !=i:
                    H[i,j] = C*2*((-1)**(i-j))/(N**2*math.sin(math.pi*(i-j)/N)*math.sin(math.pi*(i-j))/N)
    for i in range(N):
        H[i,i]=C*(1+2/N**2)/3

if N%2 == 1:
    for i in range(N):
        for j in range(N):
            if j!= i:
                H[i,j] = C*2*((-1)**(i-j)) * math.cos(math.pi*(i-j)/N)/(N**2*math.sin(math.pi*(i-j)/N)*math.sin(math.pi*(i-j)/N))
    for i in range(N):
        H[i,i] = C*(1-1/N**2)/3

for i in range(N):
    H[i,i] += V[i] 


[eig_value, eig_vector] = eigsh(H, k=n_eigs, which='SM')
print(eig_value)
#print(eig_vector)
