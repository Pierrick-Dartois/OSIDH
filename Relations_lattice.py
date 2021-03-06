import os
from Group_basis import *
from sage.all import *
from time import time

# General parameters (Colo and Kohel, OSIDH, Section 6)
n=256
t=74
l=2
d_K=-4#-3 or -4
f=l**n
d=d_K*f**2

# h_1=#Cl(O_1)
if d_K==-3:
	h_1=(l-kronecker(d_K,l))//3
else:
	h_1=(l-kronecker(d_K,l))//2

L_q=[]
q=3
for j in range(t):
	while not Mod(d_K,q).is_square() or q==l:
		q=next_prime(q)
	L_q.append(q)
	q=next_prime(q)

e=gp.qfbprimeform(d,1)
Q=[gp.qfbred(gp.qfbprimeform(d,q)) for q in L_q]

## Step 1: find a basis
t1=time()
#B,L_orders=Basis(e,Q,h,[(l,n-1)]) 
# Using this function takes too much time and is unnecessary here, see below

# Adapted to the case of class groups 
B,DL_Q_B=Basis_clg(e,Q,l,n,h_1)

t2=time()
print("Time Basis: {0} s".format(str(t2-t1)))

## Step 2: compute the matrix of discrete logarithms in this basis
L_factors=[(l,n-1)]+list(factor(h_1))
M=DL_matrix(e,Q,B,L_factors)
t3=time()
print("Time DL_matrix: {0} s".format(str(t3-t2)))

## Step 3: compute a basis of the relations lattice
if len(DL_Q_B)==1:#cyclic case
	b=[h_1*l**(n-1)]
else:
	b=[h_1*l**(n-2),l]
C=Lattice_basis(M,b)
t4=time()
print("Time Lattice_basis: {0} s".format(str(t4-t3)))
print("Total Time: {0} s".format(str(t4-t1)))

## Saving data

# Getting the coefficients of the quadratic forms appearing in B
def get_abc(q):
	r"""
	Computing the coefficients a, b, c of a quadratic form in pari type.

	INPUT: a pari/gp object gen representing a binary quadratic form (Qfb).

	OUTPUT: a tuple of python integer coefficients (a, b, c).
	"""
	L=str(q).split(",")
	a=int(L[0][4:])
	b=int(L[1])
	c=int(L[2][:-1])
	return (a,b,c)

Q_save=[get_abc(mfq) for mfq in Q]
B_save=[get_abc(mfb) for mfb in B]

# Saving the basis and the matrix
filename=os.path.join("Data_files","Relations_lattice_basis.txt")
with open(filename,"w",encoding="utf-8") as f:
	f.write("Group basis\n")
	r=len(B)
	for i in range(r):
		f.write(str(B_save[i])[1:-1]+"\n")
	f.write("Decomposition as a product of q_j:\n")
	for i in range(r):
		f.write(str(DL_Q_B[i])[1:-1]+"\n")
	f.write("Lattice basis\n")
	f.write(str(t)+"\n")
	for i in range(t):
		f.write(str(C[i])[1:-1]+"\n")


