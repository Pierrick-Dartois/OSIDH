import sys
sys.path.append("Documents/Codes/OSIDH")
from sage.all import *
#import pickle
from fpylll import *
from time import time
from OSIDH_protocol import *
from Group_basis import *

def action(chain,i,L_exp):
	r"""
	Computes the action of \prod_{j=1}^t \mfq_j^{L_exp[j]} on chain up to level i.

	INPUT:

	* chain: descending l-isogeny chain.

	* i: integer (level of computation).

	* L_exp: list of exponents of the primes.

	OUTPUT:

	i-th j-invariant of the chain (\prod_{j=1}^t \mfq_j^{L_exp[j]})*chain.
	"""

	C=Chain(chain.osidh,chain.L_j[0:i+1])
	D=C.action(L_exp)
	return D.L_j[-1]

def explore_level(L_exp,chain_pub,chain_ex,i,ind_gen,M_dl):
	r"""
	Given the right exponents L_exp such that :

	\prod_{j=1}^t \mfq_j^{L_exp[j]}*chain_pub=chain_ex up to level i

	This functions updates the list of exponent to match chain_ex up to level i+1. 
	Only works when Cl(O_n) is cyclic.

	INPUT:

	* L_exp: list of exponents.

	* chain_pub, chain_ex: descending l-isogeny chains.

	* i: level.

	* ind_gen: index of the generator of Cl(O_n) in osidh.L_mfq.

	* L_dl: matrix of the discrete logarithms of the \mfq_j in the basis [\mfq_{ind_gen}]
	(in Cl(O_n)) 

	OUTPUT:

	Returns a new chain of exponents such that:

	\prod_{j=1}^t \mfq_j^{L_exp[j]}*chain_pub=chain_ex up to level i+1
	"""

	# Test the trivial case when there is no need to change anything
	if action(chain_pub,i+1,L_exp)==chain_ex.L_j[i+1]:
		return L_exp

	# Computing h_i1=#Cl(O_{i+1}) and h_i=#Cl(O_i), which is also 
	# the exponent of the generator mfq_{ind_gen}^h_i  of Cl(O_{i+1})/Cl(O_i)
	osidh=chain_pub.osidh
	l=osidh.l
	t=osidh.t
	d_K=osidh.d_K
	if d_K==-3:
		h_i=l**(max(i-1,0))*(l-kronecker(d_K,l))//3
		h_i1=l**i*(l-kronecker(d_K,l))//3
	else:
		h_i=l**(max(i-1,0))*(l-kronecker(d_K,l))//2
		h_i1=l**i*(l-kronecker(d_K,l))//2

	# Computing the relations lattice of Cl(O_i) and LLL-reduction of the basis
	B=Lattice_basis(M_dl,[h_i1])
	B=IntegerMatrix.from_matrix(B).transpose()
	B_red=LLL.reduction(B)

	# Cardinality of Cl(O_{i+1})/Cl(O_i)
	if i==0:
		k_max=h_i1//h_i
	else:
		k_max=l

	k=1
	v_target=L_exp.copy()
	while k<k_max:
		v_target[ind_gen]+=h_i
		C=GSO.Mat(copy(B_red))
		C.update_gso()
		w=C.babai(v_target)
		w=(IntegerMatrix.from_iterable(1, t, w)*B_red)[0]
		L_exp_test=[v_target[i]-w[i] for i in range(osidh.t)]
		if action(chain_pub,i+1,L_exp_test)==chain_ex.L_j[i+1]:
			return L_exp_test
		k+=1

def attack_chain(chain_pub,chain_ex):
	r"""Computes the list of exponents e_1,..,e_t such that: 
	(\prod_{j=1}^t \mfq_j^e_j)*chain_pub=chain_ex
	Only works when Cl(O_n) is cyclic.

	INPUT:

	* chain_pub, chain_ex: descending l-isogeny chains.

	OUTPUT:

	The list of exponents [e_1,..,e_t].
	"""

	osidh=chain_pub.osidh
	l=osidh.l
	n=osidh.n
	t=osidh.t
	r=osidh.r
	d_K=osidh.d_K
	d=d_K*l**(2*n)
	L_mfq=osidh.L_mfq

	if d_K==-3:
		h_1=(l-kronecker(d_K,l))//3
		h_n=l**(n-1)*h_1
	else:
		h_1=(l-kronecker(d_K,l))//2
		h_n=l**(n-1)*h_1
	e_test=h_n//l

	# Looking for the generator of Cl(O_n)
	e=gp.qfbprimeform(d,1)
	Q=[]
	for mfq in L_mfq:
		a,b,c=get_abc(mfq)
		Q.append(gp.qfbred(gp.Qfb(a,l**n*b,l**(2*n)*c)))
	found_gen=False
	ind_gen=0
	while not found_gen:
		if Q[ind_gen]**e_test!=e:
			found_gen=True
		else:
			ind_gen+=1

	# Computing the DL of elements in Q in the basis of Q[ind_gen]
	L_factors=[(l,n-1)]+list(factor(h_1))
	L_dl=[]
	for i in range(t):
		if i==ind_gen:
			L_dl.append(1)
		else:
			L_dl.append(DL(e,[Q[ind_gen]],Q[i],L_factors,1)[0])
	M_dl=matrix(ZZ,L_dl)

	# Main loop (upgrading the exponents sucessively 
	#to match chain_ex at level i for i=1 to n)
	L_exp=[0]*t
	for i in range(n):
		L_exp=explore_level(L_exp,chain_pub,chain_ex,i,ind_gen,M_dl)
	return L_exp
		








