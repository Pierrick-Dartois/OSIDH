import sys
sys.path.append("Documents/Codes/OSIDH")
from sage.all import *
#import pickle
from fpylll import *
from time import time
from OSIDH_protocol import *
from Group_basis import *

### Part 1 - recovering the exponents when the descending l-isogeny chains are known

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

	* M_dl: matrix of the discrete logarithms of the \mfq_j in the basis [\mfq_{ind_gen}]
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


## Part 2 - recovering the descending l-isogeny chains with horizontal chains at level n

def path_to_endo(L_path1,L_path2):
	r"""
	Given two lists of horizontal chains L_path1 and L_path2 giving two path:
	phi1: E-->F and phi2: E-->F
	this functions recover the chain of prime isogenies whose composite is the endomorphism:
	hat{phi2}*phi1: E-->E

	INPUT:

	* L_path1,L_path2: list of horizontal chains.

	OUTPUT: list of isogenies L_iso with L_iso[i].codomain()==L_iso[i+1].domain() coding hat{phi2}*phi1.
	"""

	osidh=L_path1[0].osidh
	p=osidh.p
	a=osidh.F.gen()
	L_q=osidh.L_q
	N=osidh.N
	L_v=osidh.L_v

	# Choice of the first elliptic curve E
	E=EllipticCurve(j=L_path1[0].j_center)
	L_E=[E]#List of domains
	L_iso=[]#List of isogenies

	# Path 1
	for chain_hor in L_path1:
		q=L_q[chain_hor.ind_q]
		v=L_v[chain_hor.ind_q+1]
		L_j=[chain_hor.j_center]+chain_hor.L_plus+chain_hor.L_minus
		for k in range(len(L_j)-1):
			P,Q=Torsion_basis(a,p,L_E[-1],q,N,v)
			phi=L_E[-1].isogeny(P)
			if phi.codomain().j_invariant()!=L_j[k+1]:
				T=Q
				m=0
				phi=L_E[-1].isogeny(T)
				while phi.codomain().j_invariant()!=L_j[k+1] and m<=q-1:
					T+=P
					m+=1
					phi=L_E[-1].isogeny(T)
			L_iso.append(phi)
			L_E.append(phi.codomain())

	# Path 2
	t=len(L_path2)
	for j in range(t):
		chain_hor=L_chain2[t-1-j]
		q=L_q[chain_hor.ind_q]
		v=L_v[chain_hor.ind_q+1]
		L_j=[chain_hor.j_center]+chain_hor.L_plus+chain_hor.L_minus
		L_j.reverse()
		for k in range(len(L_j)-1):
			P,Q=Torsion_basis(a,p,L_E[-1],q,N,v)
			phi=L_E[-1].isogeny(P)
			if phi.codomain().j_invariant()!=L_j[k+1]:
				T=Q
				m=0
				phi=L_E[-1].isogeny(T)
				while phi.codomain().j_invariant()!=L_j[k+1] and m<=q-1:
					T+=P
					m+=1
					phi=L_E[-1].isogeny(T)
			L_iso.append(phi)
			L_E.append(phi.codomain())
	
	# Make sure we have an endomorphism
	iso=L_iso[-1].codomain().isomorphism_to(E)
	L_iso[-1].set_post_isomorphism(iso)

	return L_iso




def go_up_chain(L_chains_hor,M_dl,i,k_BKZ=4):
	r"""
	Given the data of the horizontal chains at level i:
	mfq_j^-r*Ei-->...-->Ei-->...-->mfq_j^r*Ei
	for all j in range(t), this function outputs the chains:
	mfq_j^-r*E{i-1}-->...-->E{i-1}-->...-->mfq_j^r*E{i-1}
	at level i-1.

	INPUT:

	* L_chains_hor: List of horizontal chains.

	* M_dl: matrix of the discrete logarithms of the \mfq_j in the basis [\mfq_{ind_gen}]
	(in Cl(O_n)).

	* i: level (integer >=1).

	* k_BKZ: BKZ parameter to find a short vector (integer).

	OUTPUT: List of horizontal chains at level n-1.
	"""

	# Parameters
	osidh=L_chains_hor[0].osidh
	l=osidh.l
	n=osidh.n
	t=osidh.t
	r=osidh.r
	d_K=osidh.d_K
	p=osidh.p
	a=osidh.F.gen()
	N=osidh.N
	L_v=osidh.L_v

	if d_K==-3:
		h_1=(l-kronecker(d_K,l))//3
		h_i=l**(i-1)*h_1
	else:
		h_1=(l-kronecker(d_K,l))//2
		h_i=l**(i-1)*h_1


	# Finding a short vector in the relations lattice
	B=Lattice_basis(M_dl,[h_i])
	B=IntegerMatrix.from_matrix(B).transpose()
	B_red=BKZ.reduction(B, BKZ.Param(k_BKZ))

	L_max=[]
	for k in range(t):
		L_max.append(max(max(B_red[k]),-min(B_red[k])))
	N_inf=min(L_max)
	k0=L_max.index(N_inf)
	u=M_red[k0]

	# Dividing u into u=v1-v2 with \|v_i\|_\inf<=r
	v1=[]
	v2=[]
	for k in range(t):
		if u[k]>=0:
			v1.append(r*(u[k]//r))
			v2.append(-u[k]%r)
		else:
			v1.append(-r*((-u[k])//r))
			v2.append((-u[k])%r)

	# Computing the paths Ei-->E associated to v1 and v2 
	# whose composite give the disired endomorphism
	L_path1=Action_hor_path(osidh,L_chains_hor,v1)
	L_path2=Action_hor_path(osidh,L_chains_hor,v2)

	# Computing the list of isogenies giving the desired endomorphism iota_i(beta):E-->E
	L_iso=path_to_endo(L_path1,L_path2)
	Ei=L_iso[0].domain()

	# Basis of Ei[l]
	P,Q=Torsion_basis(a,p,Ei,l,N,v)
	
	# Computing 
	R,S=P,Q
	for iso in L_iso:
		R=iso(R)
		S=iso(S)















