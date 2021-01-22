import numpy as np
from scipy import sparse
import scipy
from scipy.sparse.linalg import eigsh
from scipy.sparse.linalg import LinearOperator
import math

#parameters of variables
#number of grid points
N = 201
#mass
m = 1.0
xmax = 10.0
xmin = -10.0
L = (xmax - xmin)*(N+1)/N
x = np.linspace(xmin, xmax, N)
#for N even
C = (0.5*(math.pi*N/L)**2)/m
#number of eigenvalues
n_eigs = 5
#1D Kinetic Energy Operator:
T = np.zeros(shape=(N,N))

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
#potential function

V = np.linspace(xmin,xmax,N)
for i in range(N):
    V[i] = 0.5*V[i]**2

H = np.zeros(shape=(N,N))
H += T
for i in range(N):
    H[i,i] = H[i,i] + V[i] 

def mv(v):
    return np.dot(H,v)

A = LinearOperator((N,N),matvec=mv,dtype=float)
print(A)
w,v = scipy.sparse.linalg.eigsh(A,k=10,which="SM")
print(w)