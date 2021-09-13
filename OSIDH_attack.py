import sys
sys.path.append("Documents/Codes/OSIDH")
from sage.all import *
#import pickle
from fpylll import *
from time import time
from sage.modules.free_module_integer import IntegerLattice
from OSIDH_protocol import *
from Group_basis import *

### Part 2 - recovering the exponents when the descending l-isogeny chains are known


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

def explore_level(L_exp,chain_pub,chain_ex,i,DL_Q_B,M_dl):
	r"""
	Given the right exponents L_exp such that :

	\prod_{j=1}^t \mfq_j^{L_exp[j]}*chain_pub=chain_ex up to level i

	This function updates the list of exponent to match chain_ex up to level i+1.

	INPUT:

	* L_exp: list of exponents.

	* chain_pub, chain_ex: descending l-isogeny chains.

	* i: level.

	* DL_Q_B: "the DL of elements of the basis B_group in L_mfq (generators mfq_j)", 
	meaning that B_group[i]=\prod_{j=1}^t L_mfq[j]**DL_Q_B[i][j].  

	* M_dl: matrix of the discrete logarithms of the \mfq_j in the basis B_group
	(in Cl(O_n)).

	OUTPUT:

	Returns a new chain of exponents such that:

	\prod_{j=1}^t \mfq_j^{L_exp[j]}*chain_pub=chain_ex up to level i+1
	"""

	# Test the trivial case when there is no need to change anything
	if action(chain_pub,i+1,L_exp)==chain_ex.L_j[i+1]:
		return L_exp

	# Computing h_i1=#Cl(O_{i+1}) and h_i=#Cl(O_i)
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

	if len(DL_Q_B)==1:#Cyclic case
		# Computing the relations lattice of Cl(O_{i+1}) and LLL-reduction of the basis
		B=Lattice_basis(M_dl,[h_i1])
		B=IntegerMatrix.from_matrix(B).transpose()
		B_red=LLL.reduction(B)

		# Cardinality of ker(Cl(O_{i+1})-->>Cl(O_i))
		if i==0:
			k_max=h_i1
		else:
			k_max=l

		# Power of the generator
		if i==0:
			power=1
		else:
			power=h_i

		k=1
		v_target=L_exp.copy()
		while k<k_max:
			for j in range(t):
				v_target[j]+=power*DL_Q_B[0][j]# B_group[0]**power generates ker(Cl(O_{i+1})-->>Cl(O_i))
			C=GSO.Mat(copy(B_red))
			C.update_gso()
			w=C.babai(v_target)
			w=(IntegerMatrix.from_iterable(1, t, w)*B_red)[0]
			L_exp_test=[v_target[j]-w[j] for j in range(t)]
			if action(chain_pub,i+1,L_exp_test)==chain_ex.L_j[i+1]:
				return L_exp_test
			k+=1
	elif i>=3:#Non-cyclic case i>=3
		# Computing the relations lattice of Cl(O_{i+1}) and LLL-reduction of the basis
		B=Lattice_basis(M_dl,[h_i,l])
		B=IntegerMatrix.from_matrix(B).transpose()
		B_red=LLL.reduction(B)

		# Cardinality of ker(Cl(O_{i+1})-->>Cl(O_i))
		k_max=l
		power=h_i//l

		k=1
		v_target=L_exp.copy()
		while k<k_max:
			for j in range(t):
				v_target[j]+=power*DL_Q_B[0][j]# B_group[0]**power generates ker(Cl(O_{i+1})-->>Cl(O_i))
			C=GSO.Mat(copy(B_red))
			C.update_gso()
			w=C.babai(v_target)
			w=(IntegerMatrix.from_iterable(1, t, w)*B_red)[0]
			L_exp_test=[v_target[j]-w[j] for j in range(t)]
			if action(chain_pub,i+1,L_exp_test)==chain_ex.L_j[i+1]:
				return L_exp_test
			k+=1
	elif i==2:#Non-cyclic case i=2
		# Computing the relations lattice of Cl(O_3) and LLL-reduction of the basis
		B=Lattice_basis(M_dl,[h_i,l])
		B=IntegerMatrix.from_matrix(B).transpose()
		B_red=LLL.reduction(B)

		# Cardinality of ker(Cl(O_3)-->>Cl(O_2))
		k_max=l

		# DL of the generator of ker(Cl(O_3)-->>Cl(O_2)) in L_mfq
		DL_Q_gen=ker32(osidh,DL_Q_B)

		k=1
		v_target=L_exp.copy()
		while k<k_max:
			for j in range(t):
				v_target[j]+=DL_Q_gen[j]
			C=GSO.Mat(copy(B_red))
			C.update_gso()
			w=C.babai(v_target)
			w=(IntegerMatrix.from_iterable(1, t, w)*B_red)[0]
			L_exp_test=[v_target[j]-w[j] for j in range(t)]
			if action(chain_pub,i+1,L_exp_test)==chain_ex.L_j[i+1]:
				return L_exp_test
			k+=1
	else:#Non-cyclic case i=0,1
		DL_Q_geni,M_dli=DL_matrix_small(osidh,i)

		# k_max=#ker(Cl(O_{i+1})-->>Cl(O_i))
		# ker(Cl(O_{i+1})-->>Cl(O_i))=<gen^power>
		if i==1:
			k_max=l
			power=h_i
		else:
			k_max=h_i1
			power=1

		k=1
		while k<k_max:
			L_exp_test=[power*DL_Q_geni[j] for j in range(t)]
			if action(chain_pub,i+1,L_exp_test)==chain_ex.L_j[i+1]:
				return L_exp_test
			k+=1

def ker32(osidh,DL_Q_B):
	r"""
	Computes a generator of ker(Cl(O_3)-->>Cl(O_2)) in the non-cyclic case.

	INPUT:

	* osidh: OSIDH instanciation.

	* DL_Q_B: "the DL of elements of the basis B_group of Cl(O_3) in L_mfq (generators mfq_j)", 
	meaning that B_group[i]=\prod_{j=1}^t L_mfq[j]**DL_Q_B[i][j]. 

	OUTPUT:

	* DL_Q_gen: DL of a generator of ker(Cl(O_3)-->>Cl(O_2)) in L_mfq (generators mfq_j).
	"""

	l=osidh.l
	t=osidh.t
	d_K=osidh.d_K
	d2=d_K*l**4
	L_mfq=osidh.L_mfq

	# Computing h_1=#Cl(O_1)
	if d_K==-3:
		h_1=(l-kronecker(d_K,l))//3
	else:
		h_1=(l-kronecker(d_K,l))//2

	# Computing the image of the basis in Cl(O_2)
	e=gp.qfbprimeform(d2,1)
	g0=e
	g1=e
	for j in range(t):
		if DL_Q_B[0][j]!=0 or DL_Q_B[1][j]!=0:
			a,b,c=get_abc(L_mfq[j])
			mfq2=gp.qfbred(gp.Qfb(a,l**2*b,l**(2*2)*c))
			g0*=mfq2**DL_Q_B[0][j]
			g1*=mfq2**DL_Q_B[1][j]

	# We have ker(Cl(O_3)-->>Cl(O_2))=<g0^{ah_1}g1^b>, so we compute a and b in {0...l-1}
	h0=g0**h_1
	if h0==e:# ker=<h0>
		return [h_1*x for x in DL_Q_B[0]]
	elif g1==e:# ker=<g1>
		return DL_Q_B[1]
	else:
		a=1
		while h0**a!=g1 and a<l:
			a+=1
		DL_Q_gen=[0]*t
		for j in range(t):
			DL_Q_gen[j]=a*h_1*DL_Q_B[0][j]-DL_Q_B[1][j]
		return DL_Q_gen
	

def DL_matrix_class(osidh):
	r"""
	Computes a basis of Cl(O_n) in osidh.L_mfq and the DL 
	matrix of the primes of osidh.L_mfq in this basis.

	INPUT: 

	* osidh: an oject of the class OSIDH.

	OUTPUT: 

	* DL_Q_B: "the DL of elements of the basis B_group in L_mfq (generators mfq_j)", 
	meaning that B_group[i]=\prod_{j=1}^t L_mfq[j]**DL_Q_B[i][j].  

	* M_dl: a (1 or 2)*t integer matrix of the DL of L_mfq in B_group.
	"""

	l=osidh.l
	n=osidh.n
	t=osidh.t
	d_K=osidh.d_K
	d=d_K*l**(2*n)
	L_mfq=osidh.L_mfq

	if d_K==-3:
		h_1=(l-kronecker(d_K,l))//3
	else:
		h_1=(l-kronecker(d_K,l))//2

	# List of prime ideals (represented as quadratic forms) in O_n
	e=gp.qfbprimeform(d,1)
	Q=[]
	for mfq in L_mfq:
		a,b,c=get_abc(mfq)
		Q.append(gp.qfbred(gp.Qfb(a,l**n*b,l**(2*n)*c)))
	
	B_group,DL_Q_B=Basis_clg(e,Q,l,n,h_1)
	L_factors=[(l,n-1)]+list(factor(h_1))
	M_dl=DL_matrix(e,Q,B_group,L_factors)
	
	return DL_Q_B,M_dl

def DL_matrix_small(osidh,i):
	r"""DL_matrix_class for Cl(O_1) and Cl(O_2) (used when Cl(O_n) is not cyclic). 
	In that case, the basis of Cl(O_n) does not map to a basis of Cl(O_{1, 2}).
	Note that Cl(O_i) is always cyclic for i=1,2.

	INPUT: 

	* osidh: an oject of the class OSIDH.

	* i: integer in {1, 2}.

	OUTPUT: 

	* DL_Q_gen: "the DL of elements of the generator gen of Cl(O_i) in L_mfq", 
	meaning that gen=\prod_{j=1}^t L_mfq[j]**DL_Q_gen[j].  

	* M_dl: a 1*t integer matrix of the DL of L_mfq in the basis [gen].
	"""

	l=osidh.l
	n=osidh.n
	t=osidh.t
	d_K=osidh.d_K
	d=d_K*l**(2*i)
	L_mfq=osidh.L_mfq

	if d_K==-3:
		h_1=(l-kronecker(d_K,l))//3
	else:
		h_1=(l-kronecker(d_K,l))//2

	if h_1==1 and i==1:
		return [0]*t,matrix(ZZ,[0]*t)
	else:
		# List of prime ideals (represented as quadratic forms) in O_i
		e=gp.qfbprimeform(d,1)
		Q=[]
		for mfq in L_mfq:
			a,b,c=get_abc(mfq)
			Q.append(gp.qfbred(gp.Qfb(a,l**i*b,l**(2*i)*c)))
		
		if i==1:
			L_factors=list(factor(h_1))
			N=h_1
		else:
			L_factors=[(l,1)]+list(factor(h_1))
			N=h_1*l

		L_Nk=[N//x[0]**x[1] for x in L_factors]

		# Finding the generator: there is only one because Cl(O_i) is cyclic
		s=len(L_Nk)
		g=e
		DL_Q_gen=[0]*t
		for k in range(s):
			for j in range(t):
				h=Q[j]**L_Nk[k]
				if h**(L_factors[k][0]**(L_factors[k][1]-1))!=e:
					g*=h
					DL_Q_gen[j]+=L_Nk[k]
					break

		M_dl=DL_matrix(e,Q,[g],L_factors)

		return DL_Q_gen,M_dl


def attack_chain(chain_pub,chain_ex,DL_Q_B="omitted",M_dl="omitted"):
	r"""Computes the list of exponents e_1,..,e_t such that: 
	(\prod_{j=1}^t \mfq_j^e_j)*chain_pub=chain_ex
	Only works when Cl(O_n) is cyclic.

	INPUT:

	* chain_pub, chain_ex: descending l-isogeny chains.

	* DL_Q_B: the DL of elements of the basis B_group in L_mfq (generatos mfq_j)", 
	meaning that B_group[i]=\prod_{j=1}^t L_mfq[j]**DL_Q_B[i][j].  

	* M_dl: if not omitted, the DL matrix of the \mfq_j in 
	the basis B_group.

	OUTPUT:

	The list of exponents [e_1,..,e_t].
	"""

	osidh=chain_pub.osidh
	n=osidh.n
	t=osidh.t
	
	if DL_Q_B=="omitted":
		DL_Q_B,M_dl=DL_matrix_class(osidh)

	# Main loop (upgrading the exponents sucessively 
	#to match chain_ex at level i for i=1 to n)
	L_exp=[0]*t
	for i in range(n):
		L_exp=explore_level(L_exp,chain_pub,chain_ex,i,DL_Q_B,M_dl)
	return L_exp


### Part 1 - recovering the descending l-isogeny chains with horizontal chains at level n

def path_to_endo(L_path1,L_path2):
	r"""
	Given two lists of horizontal chains L_path1 and L_path2 giving two path:
	phi1: E-->F and phi2: E-->F
	this function recovers the chain of prime isogenies whose composite is the endomorphism:
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
	if E.cardinality()!=N:
		E=E.quadratic_twist()
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
		chain_hor=L_path2[t-1-j]
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

def generator_ideal(osidh,L_exp,i):
	r"""
	Computes the coefficients a and b of beta=a+l**i*b*theta\in O_i such that:
	beta*O_i=\prod_{j=1}^t(\mfq_j\cap O_i)^{L_exp[j]}

	INPUT:

	* osidh: an object of the class OSIDH.

	* L_exp: list of exponents (integres).

	* i: level (integer).

	OUTPUT: returns a and b (integers).
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

	return ZZ(beta[0]),ZZ(beta[1])//(l**i)

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


def go_up_chain(L_chains_hor,u,i):
	r"""
	Given the data of the horizontal chains at level i:
	mfq_j^-r*Ei-->...-->Ei-->...-->mfq_j^r*Ei
	for all j in range(t), this function outputs the chains:
	mfq_j^-r*E{i-1}-->...-->E{i-1}-->...-->mfq_j^r*E{i-1}
	at level i-1.

	INPUT:

	* L_chains_hor: List of horizontal chains.

	* u: short vector in the relations lattice of Cl(O_i) 
	(of infinity norm <=2*r)

	* i: level (integer >=1).

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

	# Dividing u into u=v1-v2 with \|v_i\|_\inf<=r
	v1=[]
	v2=[]
	for k in range(t):
		if u[k] in [-2*r,2*r]:
			v1.append(u[k]//2)
			v2.append(-u[k]//2)
		elif u[k]>=0:
			v1.append(r*(u[k]//r))
			v2.append(-(u[k]%r))
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

	# Computing a in the decomposition beta=ai+bi*l**i*theta
	ai,bi=generator_ideal(osidh,u,i)

	# Basis of Ei[l]
	P,Q=Torsion_basis(a,p,Ei,l,N,L_v[0])
	
	# Computing iota_i(bi*l**i*theta)(P,Q)
	R,S=P,Q
	for iso in L_iso:
		R=iso(R)
		S=iso(S)
	R=R-ai*P
	S=S-ai*Q

	# Finding ker(iota_i(b*l**i*theta))\cap Ei[l](=ker(T) here)
	if R.is_zero():
		T=P
	else:
		T=Q
		U=S
		k=0
		while not U.is_zero() and k<=l-1:
			U+=R
			T+=P
			k+=1

		if not U.is_zero():#sign error iso=iota_i(-ai-bi*l**i*theta)
			R+=2*ai*P
			S+=2*ai*Q
			if R.is_zero():
				T=P
			else:
				T=Q
				U=S
				k=0
				while not U.is_zero() and k<=l-1:
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
	l=osidh.l
	n=osidh.n
	t=osidh.t
	r=osidh.r
	d_K=osidh.d_K

	if d_K==-3:
		h_1=(l-kronecker(d_K,l))//3
	else:
		h_1=(l-kronecker(d_K,l))//2

	# Computing the relations lattices Li at each level i
	# and a list of short vectors vi of Li such that vi not in L{i+1}
	DL_Q_B,M_dl=DL_matrix_class(osidh)
	if len(DL_Q_B)==1:
		cyclic=True
	else:
		cyclic=False

	L_u=[]
	for i in range(1,n+1):
		h_i=h_1*l**(i-1)
		if cyclic:
			B=Lattice_basis(M_dl,[h_i])
		else:
			if i==1:
				DL_Q_gen1,M_dl1=DL_matrix_small(osidh,1)
				B=Lattice_basis(M_dl1,[h_i])
			elif i==2:
				DL_Q_gen2,M_dl2=DL_matrix_small(osidh,2)
				B=Lattice_basis(M_dl2,[h_i])
			else:
				B=Lattice_basis(M_dl,[h_i//l,l])
		B=IntegerMatrix.from_matrix(B).transpose()

		bool_end=False
		k_BKZ=4
		while not bool_end:
			B_red=BKZ.reduction(B, BKZ.Param(k_BKZ))
			k=0
			while (not bool_end) and k<=t-1:
				u=B_red[k]
				if max(max(u),-min(u))<=2*r:
					a,b=generator_ideal(osidh,u,i)
					if b%l!=0:
						bool_end=True
				k+=1
			k_BKZ+=1
		L_u.append(u)

	L_j=[L_chains_hor[0].j_center]
	L_chains_hor_i=L_chains_hor

	for i in range(n):
		L_chains_hor_i=go_up_chain(L_chains_hor_i,L_u[n-1-i],n-i)
		L_j.append(L_chains_hor_i[0].j_center)
	L_j.reverse()

	return Chain(osidh,L_j)


### Full attack

def OSIDH_attack_exe(osidh,pub_chain):
	r"""Protocol and attack execution.

	INPUT:

	* osidh: instanciation of OSIDH (parameters).

	* pub_chain: public descending l-isogeny chain used in the protocol.

	OUTPUT: prints steps and execution times. Returns None.
	"""

	print("Protocol execution:")
	
	t1=time()
	# Alice's secret key
	print("\nAlice:")
	L_exp_A=[randint(-osidh.r,osidh.r) for j in range(osidh.t)]
	print("Alice's secret key:")
	print(L_exp_A)
	chain_A=pub_chain.action(L_exp_A)
	L_chains_hor_A=[]
	for j in range(osidh.t):
		chain=chain_A
		L_plus=[]
		for k in range(osidh.r):
			chain=chain.action_prime(osidh.L_mfq[j],j)
			L_plus.append(chain.L_j[-1])
		chain=chain_A
		L_minus=[]
		for k in range(osidh.r):
			chain=chain.action_prime(osidh.L_mfq_inv[j],j)
			L_minus.append(chain.L_j[-1])
		L_chains_hor_A.append(Chain_hor(osidh,j,chain_A.L_j[-1],L_plus,L_minus))
	print("Alice's action on public chain complete.")

	# Bob's secret key
	print("\nBob:")
	L_exp_B=[randint(-osidh.r,osidh.r) for j in range(osidh.t)]
	print("Bob's secret key:")
	print(L_exp_B)
	chain_B=pub_chain.action(L_exp_B)
	L_chains_hor_B=[]
	for j in range(osidh.t):
		chain=chain_B
		L_plus=[]
		for k in range(osidh.r):
			chain=chain.action_prime(osidh.L_mfq[j],j)
			L_plus.append(chain.L_j[-1])
		chain=chain_B
		L_minus=[]
		for k in range(osidh.r):
			chain=chain.action_prime(osidh.L_mfq_inv[j],j)
			L_minus.append(chain.L_j[-1])
		L_chains_hor_B.append(Chain_hor(osidh,j,chain_B.L_j[-1],L_plus,L_minus))	
	print("Bob's action on public chain complete.")

	# Alice's computation on Bob's chain
	j_BA=Action_hor(osidh,L_chains_hor_B,L_exp_A)
	print("\nAlice's action on Bob's data complete.")

	# Bob's computation on Alice's chain
	j_AB=Action_hor(osidh,L_chains_hor_A,L_exp_B)
	print("\nBob's action on Alice's data complete.")
	t2=time()
	
	print("\nShared j-invariants coincide: {0}".format(j_AB==j_BA))
	print("Protocol execution time: {0} s".format(t2-t1))

	print("\nAttack: part 1 - recovering Alice's and Bob's chains")

	t1=time()
	chain_A_recov=recover_chain(L_chains_hor_A)
	print("Alice's chain recovered: {0}".format(chain_A_recov.L_j==chain_A.L_j))
	t1_bis=time()
	print("Time: {0} s".format(t1_bis-t1))
	chain_B_recov=recover_chain(L_chains_hor_B)
	print("Bob's chain recovered: {0}".format(chain_B_recov.L_j==chain_B.L_j))
	t2=time()
	print("Time: {0} s".format(t2-t1_bis))
	print("Total timing part 1: {0} s".format(t2-t1))

	print("\nAttack: part 2 - recovering Alice's secret exponents")
	L_exp_A_attack=attack_chain(pub_chain,chain_A_recov)
	t3=time()

	print("Timing part 2: {0} s".format(t3-t2))

	print("\nAttack: part 3 - recovering the shared secret chain")

	chain_attacked=chain_B_recov.action(L_exp_A_attack)
	t4=time()

	print("Attack is correct: {0}".format(j_AB==chain_attacked.L_j[-1]))
	print("Timing part 3: {0} s".format(t4-t3))

	print("\nTotal attack timing : {0} s".format(t4-t1))

def OSIDH_attack_stats(osidh,pub_chain,N_sample):
	r"""Provides statistics on the execution times 
	of the protocol and the attack.

	INPUT:

	* osidh: instanciation of OSIDH (parameters).

	* pub_chain: public descending l-isogeny chain used in the protocol.

	* N_sample: size of the sample to run the tests.

	OUTPUT: samples of the execution times of each step.
	"""

	L_keys_A=[]
	L_keys_B=[]
	L_protocol=[]
	L_step1a=[]
	L_step1b=[]
	L_step2=[]
	L_step3=[]
	L_attack=[]

	for run in range(N_sample):
		print("Sample number {0}".format(run))
		## Protocol

		t1=time()
		# Alice's secret key
		L_exp_A=[randint(-osidh.r,osidh.r) for j in range(osidh.t)]
		chain_A=pub_chain.action(L_exp_A)
		L_chains_hor_A=[]
		for j in range(osidh.t):
			chain=chain_A
			L_plus=[]
			for k in range(osidh.r):
				chain=chain.action_prime(osidh.L_mfq[j],j)
				L_plus.append(chain.L_j[-1])
			chain=chain_A
			L_minus=[]
			for k in range(osidh.r):
				chain=chain.action_prime(osidh.L_mfq_inv[j],j)
				L_minus.append(chain.L_j[-1])
			L_chains_hor_A.append(Chain_hor(osidh,j,chain_A.L_j[-1],L_plus,L_minus))

		# Bob's secret key
		L_exp_B=[randint(-osidh.r,osidh.r) for j in range(osidh.t)]
		chain_B=pub_chain.action(L_exp_B)
		L_chains_hor_B=[]
		for j in range(osidh.t):
			chain=chain_B
			L_plus=[]
			for k in range(osidh.r):
				chain=chain.action_prime(osidh.L_mfq[j],j)
				L_plus.append(chain.L_j[-1])
			chain=chain_B
			L_minus=[]
			for k in range(osidh.r):
				chain=chain.action_prime(osidh.L_mfq_inv[j],j)
				L_minus.append(chain.L_j[-1])
			L_chains_hor_B.append(Chain_hor(osidh,j,chain_B.L_j[-1],L_plus,L_minus))	

		# Alice's computation on Bob's chain
		j_BA=Action_hor(osidh,L_chains_hor_B,L_exp_A)

		# Bob's computation on Alice's chain
		j_AB=Action_hor(osidh,L_chains_hor_A,L_exp_B)
		t2=time()

		L_protocol.append(t2-t1)
		L_keys_A.append(L_exp_A)
		L_keys_B.append(L_exp_B)

		## Attack
		t1=time()
		chain_A_recov=recover_chain(L_chains_hor_A)
		t1_bis=time()
		chain_B_recov=recover_chain(L_chains_hor_B)
		t2=time()

		L_exp_A_attack=attack_chain(pub_chain,chain_A_recov)
		t3=time()

		chain_attacked=chain_B_recov.action(L_exp_A_attack)
		t4=time()
		
		L_step1a.append(t1_bis-t1)
		L_step1b.append(t2-t1_bis)
		L_step2.append(t3-t2)
		L_step3.append(t4-t3)
		L_attack.append(t4-t1)

	return L_keys_A,L_keys_B,L_protocol,L_step1a,L_step1b,L_step2,L_step3,L_attack
















