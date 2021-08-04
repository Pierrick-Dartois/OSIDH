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
#L_exp=[randint(-osidh.r,osidh.r) for j in range(osidh.t)]
#chain_ex=chain_pub.action(L_exp)
#print("Beginning attack: part 1 knowing the chain")
#t1=time()
#L_exp_attack=attack_chain(chain_pub,chain_ex)
#t2=time()
#print("Attack execution time: {0} s".format(t2-t1))
#print("Correctness: {0}".format(chain_ex.L_j==chain_pub.action(L_exp_attack).L_j))

L_chains_hor=[]
for j in range(t):
	chain=chain_pub
	L_plus=[]
	for k in range(r):
		chain=chain.action_prime(osidh.L_mfq[j],j)
		L_plus.append(chain.L_j[-1])
	chain=chain_pub
	L_minus=[]
	for k in range(osidh.r):
		chain=chain.action_prime(osidh.L_mfq_inv[j],j)
		L_minus.append(chain.L_j[-1])
	L_chains_hor.append(Chain_hor(osidh,j,chain_pub.L_j[-1],L_plus,L_minus))

print("\nBeginning attack: part 2 recovering the chain")
t1=time()
chain_recov=recover_chain(L_chains_hor)
t2=time()
print("Attack execution time: {0} s".format(t2-t1))
print("Correctness: {0}".format(chain_pub.L_j==chain_recov.L_j))

#ind_gen,M_dl=DL_vector_class(osidh)
#h_n=l**(n-1)
#B=Lattice_basis(M_dl,[h_n])
#B=IntegerMatrix.from_matrix(B).transpose()
#B_red=BKZ.reduction(B, BKZ.Param(4))

#L_max=[]
#for k in range(t):
	#L_max.append(max(max(B_red[k]),-min(B_red[k])))
#N_inf=min(L_max)
#k0=L_max.index(N_inf)
#u=B_red[k0]

#v1=[]
#v2=[]
#for k in range(t):
	#if u[k] in [-2*r,2*r]:
		#v1.append(u[k]//2)
		#v2.append(-u[k]//2)
	#elif u[k]>=0:
		#v1.append(r*(u[k]//r))
		#v2.append(-(u[k]%r))
	#else:
		#v1.append(-r*((-u[k])//r))
		#v2.append((-u[k])%r)

#L_path1=Action_hor_path(osidh,L_chains_hor,v1)
#L_path2=Action_hor_path(osidh,L_chains_hor,v2)
