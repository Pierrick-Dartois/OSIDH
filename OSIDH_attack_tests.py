import sys
sys.path.append("Documents/Codes/OSIDH")
from sage.all import *
from time import time
from OSIDH_protocol import *
from OSIDH_attack import *

# Parameters
n=28
t=10
l=2
r=3
d_K=-4
d=d_K*l**(2*n)

osidh=OSIDH(n,t,l,r,d_K)
chain_pub=Chain(osidh)
L_exp=[randint(-osidh.r,osidh.r) for j in range(osidh.t)]
chain_ex=chain_pub.action(L_exp)
print("Beginning attack")
t1=time()
L_exp_attack=attack_chain(chain_pub,chain_ex)
t2=time()
print("Attack execution time: {0} s".format(t2-t1))
print("Correctness: {0}".format(chain_ex.L_j==chain_pub.action(L_exp_attack).L_j))


