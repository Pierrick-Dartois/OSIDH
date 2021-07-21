import sys
sys.path.append("Documents/Codes/OSIDH")
from sage.all import *
#import pickle
from fpylll import *
from time import time
from OSIDH_protocol import *

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

def explore_level(L_exp,chain_pub,chain_ex,i,ind_gen,L_dl):
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

	* L_dl: list of the discrete logarithms of the \mfq_j in the basis of \mfq_{ind_gen} 

	OUTPUT:

	Returns a new chain of exponents such that:

	\prod_{j=1}^t \mfq_j^{L_exp[j]}*chain_pub=chain_ex up to level i+1
	"""

	# Test the trivial case when there is no need to change anything
	if action(chain_pub,i+1,L_exp)=chain_ex.L_j[i+1]:
		return L_exp

	# Computing h_i1=#Cl(O_{i+1}) and h_i=#Cl(O_i), which is also 
	# the exponent of the generator mfq_{ind_gen}^h_i  of Cl(O_{i+1})/Cl(O_i)
	osidh=chain_pub.osidh
	l=osidh.l
	d_K=osidh.d_K
	if d_K==-3:
		h_i=l**(max(i-1,0))//3*(l-kronecker(d_K,l))
		h_i1=l**i//3*(l-kronecker(d_K,l))
	else:
		h_i=l**(max(i-1,0))//2*(l-kronecker(d_K,l))
		h_i1=l**i//3*(l-kronecker(d_K,l))
	d=l**(2*i)*d_K
	
	e=gp.qfbprimeform(d,1)
	Q=[gp.qfbprimeform(d,q) for q in osidh.L_q]
	M=DL_matrix(e,Q,B,L_orders)





	# Cardinality of Cl(O_{i+1})/Cl(O_i)
	if i==0:
		k_max=h_i1//h_i
	else:
		k_max=l

	k=1
	while k<k_max:
		k+=1



