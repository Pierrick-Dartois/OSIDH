import sys
sys.path.append("Documents/Codes/OSIDH")
from sage.all import *
from OSIDH_protocol import *
from time import time

# Usual toy Parameters
n=28
t=10
l=2
r=3
# so that h\simeq l**n\simeq (2r+1)**t
d_K=-4

osidh=OSIDH(n,t,l,r,d_K)
L_q=osidh.L_q
L_mfq=osidh.L_mfq
L_mfq_inv=osidh.L_mfq_inv
C=Chain(osidh)
mfq0=L_mfq[0]
mfq1=L_mfq[1]
C=Chain(osidh)
D=Chain(osidh,C.L_j.copy())
for i in range(5):	
	t1=time()
	D=D.action_prime(mfq0,0)
	t2=time()
	print(t2-t1)
for i in range(5):	
	t1=time()
	D=D.action_prime(mfq1,1)
	t2=time()
	print(t2-t1)

E=Chain(osidh,C.L_j.copy())
for i in range(5):	
	t1=time()
	E=E.action_prime(mfq1,1)
	t2=time()
	print(t2-t1)
for i in range(5):	
	t1=time()
	E=E.action_prime(mfq0,0)
	t2=time()
	print(t2-t1)
print(D.L_j==E.L_j)
