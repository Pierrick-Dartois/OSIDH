import os, sys
#sys.path.append('Documents/Codes/OSIDH')
from sage.all import *
from time import time

dirpath="Documents/Codes/OSIDH"

# Parameters
# General parameters (OSIDH, p. 29)
n=256
t=74
l=2
d_K=-4
f=l**n
d=d_K*f**2
h=l**(n-1)# Class number
# Chosen as the first prime number p>l**(2*n)*L_q[-1]*abs(d_K) 
# and congruent to 3 mod 4 (so that y^2=x^3+x is supersingular)
p=47249115145117712178898864093677401753237285151767594263097830527675496441979179546249805026740167679179672268249176843208628683029299712491232179097440354567

# List of prime numbers
L_q=[]
q=3
for j in range(t):
	while not Mod(d_K,q).is_square() or q==l:
		q=next_prime(q)
	L_q.append(q)
	q=next_prime(q)

# We also need Phi_l
L_q=[l]+L_q

t1=time()
for i in range(t+1):
	L_indices=[]
	L_coeffs=[]
	with open(dirpath+"/Modular_polynomials/phi_j_{0}.txt".format(str(L_q[i])),"r",encoding="utf-8") as f:
		for row in f:
			L_row=row.split(" ")
			L_indices.append(L_row[0])
			c=int(L_row[1][0:-1])
			c=Mod(c,p)
			L_coeffs.append(c)
	m=len(L_indices)
	with open(dirpath+"/Modular_polynomials/phi_{0}_modp.txt".format(str(L_q[i])),"w",encoding="utf-8") as f:
		for i in range(m):
			f.write(L_indices[i]+"\n")
			f.write(str(L_coeffs[i])+"\n")
t2=time()
print(t2-t1)


