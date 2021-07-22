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
#chain_ex=chain_pub.action(L_exp)
print("Beginning attack")
t1=time()
#L_exp_attack=attack_chain(chain_pub,chain_ex)
t2=time()
print("Attack execution time: {0} s".format(t2-t1))

h_1=(l-kronecker(d_K,l))//2
L_factors=[(l,n-1)]+list(factor(h_1))
L_dl=[]
L_q=osidh.L_q
ind_gen=0
e=gp.qfbprimeform(d,1)
Q=[gp.qfbprimeform(d,q) for q in L_q]
for j in range(t):
	if j==ind_gen:
		L_dl.append(1)
	else:
		L_dl.append(DL(e,[Q[ind_gen]],Q[j],L_factors,1)[0])
M_dl=matrix(ZZ,L_dl)

i=6
h_i=h_1*l**(i-1)
di=d_K*l**(2*i)
B=Lattice_basis(M_dl,[h_i])

L_dli=[]
L_factorsi=[(l,i-1)]+list(factor(h_1))
ei=gp.qfbprimeform(di,1)
Qi=[gp.qfbprimeform(di,q) for q in L_q]
for j in range(t):
	if j==ind_gen:
		L_dli.append(1)
	else:
		L_dli.append(DL(ei,[Qi[ind_gen]],Qi[j],L_factorsi,1)[0])
M_dli=matrix(ZZ,L_dli)
Bi=Lattice_basis(M_dli,[h_i])

