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


def DL_vector_class(osidh):
	r"""
	Computes a generator of Cl(O_n) in osidh.L_mfq and the DL 
	matrix (actually a line vector) of the primes of osidh.L_mfq 
	in the basis of this generator. Only works when Cl(O_n) is 
	cyclic.


	INPUT: 

	* osidh: an oject of the class OSIDH.

	OUTPUT: 

	* ind_gen: index of the generator in osidh.L_mfq.

	* M_dl: a 1*t integer matrix.
	"""

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
	
	return ind_gen,M_dl


def attack_chain(chain_pub,chain_ex,ind_gen="omitted",M_dl="omitted"):
	r"""Computes the list of exponents e_1,..,e_t such that: 
	(\prod_{j=1}^t \mfq_j^e_j)*chain_pub=chain_ex
	Only works when Cl(O_n) is cyclic.

	INPUT:

	* chain_pub, chain_ex: descending l-isogeny chains.

	* ind_gen: if not omitted, index of a generator of Cl(O_n) 
	in the list of the \mfq_j.

	* M_dl: if not omitted, the DL matrix of the \mfq_j in 
	the basis [\mfq_{ind_gen}].

	OUTPUT:

	The list of exponents [e_1,..,e_t].
	"""

	osidh=chain_pub.osidh
	n=osidh.n
	t=osidh.t
	
	if ind_gen=="omitted":
		ind_gen,M_dl=DL_vector_class(osidh)

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

def integral_part(osidh,L_exp,i):
	r"""
	Computes the integral part a of beta=a+l**i*b*theta\in O_i such that:
	beta O_i=\prod_{j=1}^t(\mfq_j\cap O_i)^{L_exp[j]}

	INPUT:

	* osidh: an object of the class OSIDH.

	* L_exp: list of exponents (integres).

	* i: level (integer).

	OUTPUT: returns a (integer).
	"""

	t=osidh.t
	L_mfq=osidh.L_mfq
	L_mfq_inv=osidh.L_mfq_inv
	d_K=osidh.d_K
	l=osidh.l

	PQ=QQ['x']
	x=PQ.gen()
	L_ideals=[]
	L_ideals_inv=[]
	if d_K==-3:
		K=NumberField(x**2+x+1,"theta")
		theta=K.gen()
		for j in range(t):
			a,b,c=get_abc(L_mfq[j])
			L_ideals.append(K.ideal([a,(1-b)//2+theta]))
			a,b,c=get_abc(L_mfq_inv[j])
			L_ideals_inv.append(K.ideal([a,(1-b)//2+theta]))
	else:
		K=NumberField(x**2+1,"theta")
		theta=K.gen()
		for j in range(t):
			a,b,c=get_abc(L_mfq[j])
			L_ideals.append(K.ideal([a,-b//2+theta]))
			a,b,c=get_abc(L_mfq_inv[j])
			L_ideals_inv.append(K.ideal([a,-b//2+theta]))

	I=K.ideal(1)
	for j in range(t):
		if L_exp[j]>=0:
			I*=L_ideals[j]**L_exp[j]
		else:
			I*=L_ideals_inv[j]**(-L_exp[j])

	beta=I.gens_reduced()[0]
	while beta[1]%(l**i)!=0:
		beta*=theta

	return beta[0]

def chain_hor_up(chain_hor,j):
	r"""Given the horizontal chain at level i:
	mfq^-r*Ei-->...-->Ei-->...-->mfq^r*Ei
	and the elliptic curve E{i-1} (represented by j).
	This function computes the chain at level i-1:
	mfq^-r*E{i-1}-->...-->E{i-1}-->...-->mfq^r*E{i-1}	
	

	INPUT:

	* chain_hor: the horizontal chain:
	mfq^-r*Ei-->...-->Ei-->...-->mfq^r*Ei

	* j: the j-invariant of E{i-1}

	OUTPUT: The horizontal chain:
	mfq^-r*E{i-1}-->...-->E{i-1}-->...-->mfq^r*E{i-1}
	"""

	osidh=chain_hor.osidh
	L_phi=osidh.L_phi
	r=osidh.r
	Fz=osidh.Fz

	phi_l=L_phi[0]
	phi_q=L_phi[chain_hor.ind_q+1]
	
	L_plus=[j]
	L_minus=[j]

	k=1
	while k<=r:
		L_eval=[]
		for g in phi_l:
			L_eval.append(g(chain_hor.L_plus[k-1]))
		f_l=Fz(L_eval)
		L_eval=[]
		for g in phi_q:
			L_eval.append(g(L_plus[-1]))
		f_q=Fz(L_eval)
		f=gcd(f_l,f_q)
		L_roots=f.roots(multiplicities=False)
		L_plus.append(L_roots[0])
		k+=1
	k=1
	while k<=r:
		L_eval=[]
		for g in phi_l:
			L_eval.append(g(chain_hor.L_minus[k-1]))
		f_l=Fz(L_eval)
		L_eval=[]
		for g in phi_q:
			L_eval.append(g(L_minus[-1]))
		f_q=Fz(L_eval)
		f=gcd(f_l,f_q)
		L_roots=f.roots(multiplicities=False)
		L_minus.append(L_roots[0])
		k+=1
	return Chain_hor(osidh,chain_hor.ind_q,j,L_plus[1::],L_minus[1::])


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

	# Computing a in the decomposition beta=a+b*l**i*theta
	a=integral_part(osidh,u,i)

	# Basis of Ei[l]
	P,Q=Torsion_basis(a,p,Ei,l,N,v)
	
	# Computing iota_i(b*l**i*theta)(P,Q)
	R,S=P,Q
	for iso in L_iso:
		R=iso(R)
		S=iso(S)
	R=R-a*P
	S=S-a*Q

	# Finding ker(iota_i(b*l**i*theta))\cap Ei[l](=ker(T) here)
	if R.is_zero():
		T=P
	else:
		T=Q
		U=S
		k=0
		while not U.is_zero() and k<=q-1:
			U+=R
			T+=P
			k+=1

	# Computing E{i-1}
	phi_dual=Ei.isogeny(T)
	jim1=phi_dual.codomain().j_invariant()

	# Computing the horizontal chains at level i-1
	L_chains_hor_im1=[]
	for j in range(t):
		L_chains_hor_im1.append(chain_hor_up(L_chains_hor[j],jim1))
	return L_chains_hor_im1

def recover_chain(L_chains_hor):
	r"""
	Given the data of the horizontal chains at level n:
	mfq_j^-r*En-->...-->En-->...-->mfq_j^r*En
	This function recovers the whole descending l-isogeny chain:
	E0-->...-->En

	INPUT:

	* L_chains_hor: List of horizontal chains.
	
	OUTPUT: The associated descending l-isogeny chain.
	"""

	
	osidh=L_chains_hor[0].osidh
	n=osidh.n

	ind_gen,M_dl=DL_vector_class(osidh)

	L_j=[L_chains_hor[0].j_center]
	L_chains_hor_i=L_chains_hor

	for i in range(n-1):
		L_chains_hor_i=go_up_chain(L_chains_hor_i,M_dl,n-i)
		L_j.append(L_chains_hor_i[0].j_center)
	L_j.append(osidh.E0.j_invariant())

	return Chain(osidh,L_j)






















