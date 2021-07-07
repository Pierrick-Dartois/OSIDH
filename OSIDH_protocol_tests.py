import sys
sys.path.append("Documents/Codes/OSIDH")
from sage.all import *
from OSIDH_protocol import *
from time import time

## Parameters
# General parameters (toy parameters)
n=28
t=10
l=2
r=3
# so that h\simeq l**n\simeq (2r+1)**t
d_K=-4
f=l**n
d=d_K*f**2
h=l**(n-1)# Class number

# List of prime numbers
L_q=[]
q=3
for j in range(t):
	while not Mod(d_K,q).is_square() or q==l:
		q=next_prime(q)
	L_q.append(q)
	q=next_prime(q)


# Chosen as the first prime number p>l**(2*n)*L_q[-1]*abs(d_K) 
# and congruent to 3 mod 4 (so that y^2=x^3+x is supersingular)
p=next_prime(L_q[-1]*2**(2*n)*4)
while Mod(p,4)!=3: 
	p=next_prime(p)

F=GF(p**2,"a")
R=PolynomialRing(F,2,["x","y"])
x,y=R.gens()
S=PolynomialRing(F,"z")

## Test Phi(R,q,x,y)
print("Testing the function Phi")

# Test 1: Phi_2
try:
	f=Phi(R,2,x,y)
	print("Test 1: success")
	print("Phi_2(x,y) = {0}".format(str(f)))
except:
	print("Test 1: failure")

# Test 2: Phi_q for q prime in L_q
try:
	t1=time()
	f=Phi(R,L_q[0],x,y)
	t2=time()
	print("Test 2: success")
	print("Time of test 2: {0} s".format(t2-t1))
except:
	print("Test 2: failure")

#Test 3: not in the database
try:
	f=Phi(R,997,x,y)
	print("Test 3: failure")
except:
	print("Test 3: success")


## Test Chain(R,S,j0,l,n)
print("\nTesting the function Chain")

# Test 1:
try:
	L=Chain(R,S,1728,l,n)
	print("Test 1: sucess")
except:
	print("Test 1: failure")


