import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
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
C = (0.5*(3.14*N/L)**2)/m
#number of eigenvalues
n_eigs = 5
#potential function
#def V_function(x):
#        V=.5*x**2
#        return V
V = np.linspace(xmin,xmax,N)
for i in range(N):
    V[i] = .5*V[i]**2

#dx = x[1] - x[0]
#V = V_function(x)
    #Hamiltonian Operator matrix:
H = sparse.eye(N, N, format='lil')*0
    # implement the numerical derivative

if N%2 == 0:
    for i in range(N):
        for j in range(N):
	        if j !=i:
                    H[i,j] = C*2*((-1)**(i-j))/(N**2*math.sin(3.14*(i-j)/N)*math.sin(3.14*(i-j))/N)
    for i in range(N):
        H[i,i]=C*(1+2/N**2)/3

if N%2 == 1:
    for i in range(N):
        for j in range(N):
            if j!= i:
                H[i,j] = C*2*((-1)**(i-j)) * math.cos(3.14*(i-j)/N)/(N**2*math.sin(3.14*(i-j)/N)*math.sin(3.14*(i-j)/N))
    for i in range(N):
        H[i,i] = C*(1-1/N**2)/3

for i in range(N):
    H[i,i] += V[i] 
[eig_value, eig_vector] = eigsh(H, k=n_eigs, which='SM')
print(eig_value)
#print(eig_vector)
