import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
import math

#parameters of variables
#number of grid points
N = 201
#mass
m = 1.0
xmax = 10.0
xmin = -10.0
#L is not simply difference between xmax and xmin, because box looks like this:
# | x  x  x  x  x  x |
#Thus, L is equal to (N+1)(xmax-xmin)/N

L = (xmax - xmin)*(N+1)/N
x = np.linspace(xmin, xmax, N)
#for N even
C = (0.5*(math.pi*N/L)**2)/m
#number of eigenvalues
n_eigs = 5
#potential function
#def V_function(x):
#        V=.5*x**2
#        return V
V = np.linspace(xmin,xmax,N)
for i in range(N):
    V[i] = 0.5*V[i]**2

#dx = x[1] - x[0]
#V = V_function(x)

    #1D Kinetic Energy Operator:
T = np.zeros(shape=(N,N))
    # implement the DVR kinetic Energy Operator derivative

Tfac=math.pi*float(N)/L
Tfac=0.5*Tfac*Tfac/m
N_inv=1.0/float(N)
N_invsq=N_inv*N_inv

if N%2 == 0:
    for i in range(N):
        for j in range(N):
	        if j !=i:
                    idel = i-j
                    isgn = (-1)**idel
                    sinfac = np.sin(math.pi*float(idel)*N_inv)
                    sinfac = sinfac*sinfac
                    T[i,j] = 2.0*float(isgn)*N_invsq/sinfac
    for i in range(N):
        T[i,i]=(1.0+2.0*N_invsq)/3.0

if N%2 == 1:
    for i in range(N):
        for j in range(N):
            if j!= i:
		idel = i-j
		isgn=(-1)**idel
		sinfac = np.sin(math.pi*float(idel)*N_inv)
		cosfac = np.cos(math.pi*float(idel)*N_inv)
		sinfac = sinfac * sinfac
		T[i,j] = 2.0*float(isgn)*N_invsq*cosfac/sinfac
    for i in range(N):
        T[i,i] = (1.0-1.0*N_invsq)/3.0

T*=Tfac
#For testing - construct and store full Hamiltonian
H = np.zeros(shape=(N,N))
H += T
for i in range(N):
    H[i,i] = H[i,i] + V[i] 


#[eig_value, eig_vector] = eigsh(H, k=n_eigs, which='SM')
#print(eig_value)

w = np.linalg.eigvalsh(H)
print(w)
#print(eig_vector)
