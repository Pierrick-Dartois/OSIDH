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
C=Chain(osidh)

t1=time()
L_chains_hor=[]
for j in range(0):
	chain=C
	L_plus=[]
	for k in range(osidh.r):
		chain=chain.action_prime(osidh.L_mfq[j],j)
		L_plus.append(chain.L_j[-1])		
	chain=C
	L_minus=[]
	for k in range(osidh.r):
		chain=chain.action_prime(osidh.L_mfq_inv[j],j)
		L_minus.append(chain.L_j[-1])
	L_chains_hor.append(Chain_hor(osidh,j,chain.L_j[-1],L_plus,L_minus))
t2=time()
print("Time horizontal chains: {0} s".format(t2-t1))

t1=time()
#chain_test=L_chains_hor[0].action_step(1,L_chains_hor[1].L_plus[0],-3)
t2=time()
print("Time action_step: {0} s".format(t2-t1))

t1=time()
#chain_test=L_chains_hor[0].action_chain(L_chains_hor[1],-3,2)
t2=time()
print("Time action_chain: {0} s".format(t2-t1))

t1=time()
OSIDH_exe(osidh,C)
t2=time()
print("OSIDH execution time: {0} s".format(t2-t1))
