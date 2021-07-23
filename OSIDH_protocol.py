from sage.all import *
import pickle
from time import time

dirpath="Documents/Codes/OSIDH"

## Intermediary functions

def get_abc(q):
	r"""
	Computing the coefficients a, b, c of a qudratic form in pari type.

	INPUT: a pari/gp object gen representing a binary quadratic form (Qfb).

	OUTPUT: a tuple of sage integer coefficients (a, b, c).
	"""
	L=str(q).split(",")
	a=ZZ(L[0][4:])
	b=ZZ(L[1])
	c=ZZ(L[2][:-1])
	return (a,b,c)

def RandomFieldElt(a,p):
	r"""
	Replacement for the instruction F.random_element() where F is the 
	quadratic extension of Fp, which is very inefficient.


	INPUT:

	* a: generator of F.

	* p: characteristic of F.

	OUTPUT: a random element in F.
	"""

	return randint(0,p-1)+a*randint(0,p-1)

def RandomECElt(a,p,E):
	r"""
	Replacement for E.random_element() and E.random_point() where E is an
	elliptic curve defined over a quadractic extension F of Fp, which are very
	inefficient.

	INPUT: 

	* a: generator of F.

	* p: characteristic of F.

	* E: an elliptic curve defined over F.

	OUTPUT: a random element in E\{(0:1:0)}. The output distribution is not uniform
	because it is not necessary.
	"""

	x=RandomFieldElt(a,p)
	A=E.a4()
	B=E.a6()
	y2=x**3+A*x+B
	while not y2.is_square():
		x=RandomFieldElt(a,p)
		y2=x**3+A*x+B
	y=y2.sqrt(extend=False,all=False)
	return E((x,y))

def Torsion_basis(a,p,E,q,N,v):
	r"""
	Finds a basis of the q-torsion subgroup E[q] of the elliptic curve E, 
	provided that all its points are rational.

	INPUT:

	* a: generator of a finite field F (quadratic extension of Fp).

	* p: characteristic of F.

	* E: elliptic curve defined over F.

	* q: (prime) integer.

	* N: cardinality of E.

	OUTPUT:

	A basis of E[q].
	"""

	m=N//q**v
	if v==2:
		P=m*RandomECElt(a,p,E)
		while P.is_zero():
			P=m*RandomECElt(a,p,E)

		Q=m*RandomECElt(a,p,E)
		while P.weil_pairing(Q,q).is_one():
			Q=m*RandomECElt(a,p,E)
	else:
		P=m*RandomECElt(a,p,E)
		while P.is_zero():
			P=m*RandomECElt(a,p,E)
		# Make sure that P has order q
		P1=q*P
		while not P1.is_zero():
			P=P1
			P1=q*P1

		Q=m*RandomECElt(a,p,E)
		# Make sure that Q has order q
		Q1=q*Q
		while not Q1.is_zero():
			Q=Q1
			Q1=q*Q1
		while P.weil_pairing(Q,q).is_one():
			Q=m*RandomECElt(a,p,E)
			# Make sure that Q has order q
			Q1=q*Q
			while not Q1.is_zero():
				Q=Q1
				Q1=q*Q1
	return (P,Q)


## Base class with OSIDH parameters
class OSIDH:
	def __init__(self,n,t,l,r,d_K=-4):
		r"""
		Instanciaates the parameters of a session of OSIDH

		INPUT:

		* n: depth of the descending l-isogeny chains.

		* t: number of prime ideals q_j.

		* l: prime number (degree of the descending isogenies).

		* r: integral range for the exponents of the q_j.

		* d_K: fundamental discriminant such that h(d_K)=1 (d_K=-3 or -4).
		"""
		
		if not ZZ(l).is_prime():
			raise ValueError("l should be a prime")
		if ZZ(d_K) not in [-3,-4]:
			raise ValueError("d_K should be -3 or -4")

		self.n=n
		self.t=t
		self.l=ZZ(l)
		self.r=r
		self.d_K=ZZ(d_K)

		# List of prime numbers q_j
		L_q=[]
		q=2
		prod=l
		for j in range(t):
			while kronecker(d_K,q)==-1 or q==l:
				q=next_prime(q)
			L_q.append(q)
			prod*=q
			q=next_prime(q)
		self.L_q=L_q

		# List of prime ideals lying above the primes of L_q in O_K 
		# (represented as forms of discriminant d_K)
		self.L_mfq=[gp.qfbprimeform(d_K,q) for q in L_q]
		self.L_mfq_inv=[mfq**(-1) for mfq in self.L_mfq]

		# Make sure that p is big enough to ensure that we con compute the chains efficiently
		threshold=l**(2*n)*L_q[-1]*-d_K
		while prod<threshold:
			prod*=l


		# Prime p: p%3==2 if d_K=-3 and p%4==3 if d_K=-4 to ensure that 
		# K-oriented supersingular elliptic curves/Fp^2 exist
		# p\pm 1 divisible by l and the q_j so that torsion points are Fp^2-rational 
		if d_K==-3:
			f=1
			if (f*prod)%3==2:
				f+=1
			if (f*prod)%3==1:
				p=f*prod+1
			else:
				p=f*prod-1
			while not p.is_prime():
				f+=1
				if (f*prod)%3==2:
					f+=1
				if (f*prod)%3==1:
					p=f*prod+1
				else:
					p=f*prod-1
		else:
			f=1
			if (f*prod)%4 in [1,3]:
				f+=1
			if (f*prod)%4==2:
				p=f*prod+1
			else:
				p=f*prod-1
			while not p.is_prime():
				f+=1
				if (f*prod)%4 in [1,3]:
					f+=1
				if (f*prod)%4==2:
					p=f*prod+1
				else:
					p=f*prod-1

		self.p=p
		# Cordinality of the first elliptic curve self.E0
		self.N=(f*prod)**2
		# List of valuations
		L_v=[2]*(t+1)
		m=f*prod
		m=m//l
		while m%l==0:
			m=m//l
			L_v[0]+=2
		for i in range(t):
			m=m//L_q[i]
			while m%L_q[i]==0:
				m=m//L_q[i]
				L_v[i+1]+=2
		self.L_v=L_v


		# Field of definition and polynomial rings
		self.F=GF(p**2,"a",proof="False")
		self.Fz=PolynomialRing(self.F,"z",sparse=True)

		# Modular polynomials
		self.L_phi=[]
		for q in [l]+L_q:
			try:
				with open(dirpath+"/Modular_polynomials/phi_j_{0}.txt".format(str(q)),"r",encoding="utf-8") as file:
					L_P=[[0]*(q+2) for i in range(q+2)]
					for row in file:
						L_row=row.split(" ")
						L_ind=L_row[0][1:-1].split(",")
						u,v=int(L_ind[0]),int(L_ind[1])
						c=ZZ(L_row[1][0:-1])
						if u!=v:
							L_P[u][v]=c
							L_P[v][u]=c
						else:
							L_P[u][u]=c
					phi=[self.Fz(elt) for elt in L_P]
					self.L_phi.append(phi)
			except:
				raise ValueError("Prime q not in the database")

		# Origin curve
		a=self.F.gen()
		if d_K==-3:
			if p==f*prod-1:
				self.E0=EllipticCurve(self.F,[0,1])
			else:
				b=RandomFieldElt(a,p)
				while b.is_square():
					b=RandomFieldElt(a,p)
				self.E0=EllipticCurve(self.F,[0,b**3])
		else:
			if p==f*prod-1:
				self.E0=EllipticCurve(self.F,[1,0])
			else:
				b=RandomFieldElt(a,p)
				while b.is_square():
					b=RandomFieldElt(a,p)
				self.E0=EllipticCurve(self.F,[b**2,0])

		# Defining the orientation in E0 as follows:
		# * If d_K==-3, theta=(-1+sqrt(3)*i)/2 maps to (x,y)|-->(zeta*x,y)
		# with zeta, primitive third root of unity
		# * If d_K==-4, theta=i maps to (x,y)|-->(-x,zeta*y) with zeta**2=-1
		if d_K==-3:
			self.zeta=self.F.zeta(3)
		else:
			self.zeta=self.F(-1).sqrt(extend=False,all=False)


	def save(self,filename):
		with open(filename,"wb") as f:
			pickle.dump(self,f)


## Class of descending l-isogeny chains for a given set of parameters
class Chain:
	def __init__(self,osidh,L_j=[]):
		r"""
		Instanciates a descending l-isogeny chain of length n starting from the origin j-invariant.

		INPUT:

		* osidh: An OSIDH object (instanciates general parameters).

		* L_j: list of j-invariants.

		OUTPUT:

		If L_j is empty, instanciates a random chain. Otherwise, self.L_j=L_j.
		"""

		self.osidh=osidh
		if len(L_j)>0:
			self.L_j=L_j
		else:
			self.L_j=[]
			if osidh.d_K==-3:
				self.L_j.append(osidh.F(0))
			else:
				self.L_j.append(osidh.F(1728))

			phi_l=osidh.L_phi[0]
			Fz=osidh.Fz

			if osidh.n>=1:
				# Choice of the second j-invariant
				L_eval=[]
				for g in phi_l:
					L_eval.append(g(self.L_j[0]))
				f=Fz(L_eval)
				L_roots=f.roots(multiplicities=False)
				m=len(L_roots)
				i=randint(0,m-1)
				while L_roots[i]==self.L_j[0]:
					i=randint(0,m-1)
				self.L_j.append(L_roots[i])

				# Choice of the other j-invariants
				k=2
				while k<=osidh.n:
					L_eval=[]
					for g in phi_l:
						L_eval.append(g(self.L_j[k-1]))
					f=Fz(L_eval)
					L_roots=f.roots(multiplicities=False)
					m=len(L_roots)
					i=randint(0,m-1)
					while L_roots[i]==self.L_j[k-2]:
						i=randint(0,m-1)
					self.L_j.append(L_roots[i])
					k+=1

	def action_torsion(self,mfq,ind_q,i):
		r"""Computes the action of the prime ideal represented by the form 
		mfq at level i.
		Used to remove the ambiguity with the action by the conjugate of mfq at level i. 
		Does not use modular polynomials.
		
		INPUT:

		* mfq: prime form representing a prime ideal.
		
		* ind_q: index of mfq in self.osidh.L_mfq or self.osidh.L_mfq_inv.

		* i: level.

		OUTPUT:

		The j-invariant of mfq*Ei. 
		"""
		
		# General OSIDH parameters
		p=self.osidh.p
		N=self.osidh.N
		L_v=self.osidh.L_v
		a=self.osidh.F.gen()

		# Recovering the isogeny chain up to depth i form the j-invariants
		L_E=[self.osidh.E0]
		l=self.osidh.l
		L_iso=[]
		L_iso_dual=[]
		for j in range(i):
			# Torsion subgroup Ej[l]
			Pj,Qj=Torsion_basis(a,p,L_E[-1],l,N,L_v[0])
			# Trying all cyclic subgroups of Ej[l]
			phij=L_E[-1].isogeny(Pj)
			if phij.codomain().j_invariant()!=self.L_j[j+1]:
				T=Qj
				phij=L_E[-1].isogeny(T)
				k=0
				while phij.codomain().j_invariant()!=self.L_j[j+1] and k<=l-1:
					T+=Pj
					phij=L_E[-1].isogeny(T)
					k+=1
			
			# Computing the dual isogeny of phij (since the sage dual method is very slow)
			Ej1=phij.codomain()
			Rj,Sj=phij(Pj),phij(Qj)
			if Rj.is_zero():
				phij_dual=Ej1.isogeny(Sj)
			else:
				phij_dual=Ej1.isogeny(Rj)
			# Make sure that the codomain of phij_dual is the domain of phij
			iso=phij_dual.codomain().isomorphism_to(L_E[-1])
			phij_dual.set_post_isomorphism(iso)

			L_iso.append(phij)
			L_iso_dual.append(phij_dual)
			L_E.append(Ej1)

		# Find lamb such that the ideal associated to mfq = [q,theta-lamb]
		# where theta generates the order O_K
		q,b,c=get_abc(mfq)
		if self.osidh.d_K==-3:
			lamb=(b-1)//2
		else:
			lamb=b//2

		# Find a basis of Ei[q]
		P,Q=Torsion_basis(a,p,L_E[-1],q,N,L_v[ind_q+1])

		# Image of the basis by dual isogenies
		R,S=P,Q
		for j in range(i):
			R=L_iso_dual[i-1-j](R)
			S=L_iso_dual[i-1-j](S)

		# Image of R,S by theta-lamb
		zeta=self.osidh.zeta
		xR,yR=R.xy()
		xS,yS=S.xy()
		if self.osidh.d_K==-3:
			R1=L_E[0](zeta*xR,yR)
			S1=L_E[0](zeta*xS,yS)
		else:
			R1=L_E[0](-xR,zeta*yR)
			S1=L_E[0](-xS,zeta*yS)
		R=R1-lamb*R
		S=S1-lamb*S

		# Image of R, S by the isogeny chain E0 -->...--> Ei
		for j in range(i):
			R=L_iso[j](R)
			S=L_iso[j](S)

		# Computing a generator T of Ei[mfq]
		if R.is_zero():
			T=P
		else:
			U=R
			k=1
			while U!=S:
				U+=R
				k+=1
			T=Q-k*P

		# Computing mfq*Ei:=Ei/<T>
		phi=L_E[-1].isogeny(T)
		return phi.codomain().j_invariant()

	def action_prime(self,mfq,ind_q):
		r"""Computes the action of the prime ideal represented by the form mfq
		on self using modular polynomials. Uses action_torsion in case of 
		ambiguity (when there are two roots).

		INPUT:

		* mfq: prime form representing a prime ideal.

		* ind_q: index of an ideal mfq in self.osidh.L_mfq or self.osidh.L_mfq_inv.

		OUTPUT:

		The chain obtained via the action of mfq=self.osidh.L_mfq[ind_q] on self.
		"""

		phi_l=self.osidh.L_phi[0]
		phi_q=self.osidh.L_phi[ind_q+1]
		Fz=self.osidh.Fz

		LF_j=[self.L_j[0]]
		for i in range(len(self.L_j)-1):
			L_eval=[]
			for g in phi_l:
				L_eval.append(g(LF_j[i]))
			f_l=Fz(L_eval)
			L_eval=[]
			for g in phi_q:
				L_eval.append(g(self.L_j[i+1]))
			f_q=Fz(L_eval)
			f=gcd(f_l,f_q)
			L_roots=f.roots(multiplicities=False)
			if len(L_roots)==1:
				LF_j.append(L_roots[0])
			elif len(L_roots)>=2:
				LF_j.append(self.action_torsion(mfq,ind_q,i+1))
			else:
				raise ValueError("Modular polynomials breakdown")

		return Chain(self.osidh,LF_j)

	def action(self,L_exp):
		r"""Computes the action of an ideal of O_n written as a product of the mfq_j.

		INPUT:

		* L_exp: list of exponents in the interval [-r,r].

		OUTPUT:

		The chain obtained via the action of the ideal on self.
		"""

		out_chain=Chain(self.osidh,self.L_j.copy())# Deep copy of self
		L_mfq=self.osidh.L_mfq
		L_mfq_inv=self.osidh.L_mfq_inv

		for j in range(self.osidh.t):
			if L_exp[j]>=0:
				for k in range(L_exp[j]):
					out_chain=out_chain.action_prime(L_mfq[j],j)
			else:
				for k in range(-L_exp[j]):
					out_chain=out_chain.action_prime(L_mfq_inv[j],j)
		return out_chain


## Class of horizontal isogeny chains at the bottom for a given set of parameters
class Chain_hor:
	def __init__(self,osidh,ind_q,j_center,L_plus,L_minus):
		r"""Instanciates a horizontal q-isogeny chain at the bottom for a given prime q.
		
		INPUT:

		* osidh: An OSIDH object (instanciates general parameters).

		* ind_q: index of q in osidh.L_q.

		* j_center: central j-invariant of the chain.

		* L_plus: list of j-invariants with positive exponents.
		mfq*E-->...-->mfq^r*E

		* L_minus: list of j-invariants with negative exponents.
		mfq^-r*E-->...-->mfq^-1*E

		OUTPUT: 

		A horizontal q-isogeny chain.
		"""

		self.osidh=osidh
		self.ind_q=ind_q
		self.j_center=j_center
		self.L_plus=L_plus
		self.L_minus=L_minus

	def action_step(self,ind_q,j,e):
		r""" Computes the action of the prime ideal self.osidh.L_mfq(_inv)[ind_q]
		relating self.j_center with j on the chain self up to exponent e.

		INPUT:

		* ind_q: index of the prime in osidh.L_q.

		* j: j-invariant.

		* e: exponent.

		OUTPUT: 

		The horizontal self.osidh.L_q[ind_q]-isogeny chain obtained by pusching the action
		on self up to exponent e. 
		"""

		phi_q1=self.osidh.L_phi[self.ind_q+1]# Modular polynomial of the chain self
		phi_q2=self.osidh.L_phi[ind_q+1]# Modular polynomial of the action by self.osidh.L_mfq(_inv)[ind_q]
		Fz=self.osidh.Fz

		L_plus=[j]
		L_minus=[j]

		if e>=0:
			k=1
			while k<=e:
				L_eval=[]
				for g in phi_q1:
					L_eval.append(g(L_plus[-1]))
				f_q1=Fz(L_eval)
				L_eval=[]
				for g in phi_q2:
					L_eval.append(g(self.L_plus[k-1]))
				f_q2=Fz(L_eval)
				f=gcd(f_q1,f_q2)
				L_roots=f.roots(multiplicities=False)
				L_plus.append(L_roots[0])
				k+=1
		else:
			k=1
			while k<=-e:
				L_eval=[]
				for g in phi_q1:
					L_eval.append(g(L_minus[-1]))
				f_q1=Fz(L_eval)
				L_eval=[]
				for g in phi_q2:
					L_eval.append(g(self.L_minus[k-1]))
				f_q2=Fz(L_eval)
				f=gcd(f_q1,f_q2)
				L_roots=f.roots(multiplicities=False)
				L_minus.append(L_roots[0])
				k+=1
		return Chain_hor(self.osidh,self.ind_q,j,L_plus[1::],L_minus[1::])

	def action_chain(self,chain,e,f):
		r""" Given two chains:

		* C1: E-->...-->mfq^e*E (self)

		* C2: E-->...-->mfq'^f*E

		This method outputs mfq^e*E-->...-->mfq^e*mfq'^f*E.

		INPUT:

		* self: a chain containing C1 (with center=E).

		* chain: a chain containing C2 (with center=E).

		* e,f: integer exponents in C1 and C2.

		OUTPUT: The chain mfq^e*E-->...-->mfq^e*mfq'^f*E.
		"""

		if e==0:
			if f>=0:
				return Chain_hor(chain.osidh,chain.ind_q,chain.j_center,chain.L_plus[0:f],[])
			else:
				return Chain_hor(chain.osidh,chain.ind_q,chain.j_center,[],chain.L_minus[0:-f])
		elif e>0:
			chain_out=chain
			for k in range(e):
				chain_out=chain_out.action_step(self.ind_q,self.L_plus[k],f)
		else:
			chain_out=chain
			for k in range(-e):
				chain_out=chain_out.action_step(self.ind_q,self.L_minus[k],f)
		return chain_out


def Action_hor_path(osidh,L_chains_hor,L_exp):
	r"""Given the horizontal chains at the bottom and the exponents, this function performs the
	action of the exponents on the central j-invariant of all the chains.

	INPUT:

	* osidh: OSIDH instanciation.

	* L_chains_hor: list of horizontal chains:
	mfq_j^-r*E_n-->...-->E_n -->...-->mfq_0^e_0*E_n

	* L_exp: list of exponents (integers in [-r,r])

	OUTPUT:

	Returns the list of chains:
	E_n -->...-->mfq_0^e_0*E_n
	mfq_0^e_0*E_n -->...-->mfq_0^e_0*mfq_1^e_1*E_n
	...
	mfq_0^e_0*...*mfq_{t-1}^e_{t-1}*E_n -->...-->mfq_0^e_0*...*mfq_t^e_t*E_n
	"""
	t=osidh.t
	L_chains_path=[]
	# List of chains: at step j, this list contains the chains:
	# E_n -->...-->mfq_0^e_0*E_n
	# mfq_0^e_0*E_n -->...-->mfq_0^e_0*mfq_1^e_1*E_n
	# ...
	# mfq_0^e_0*...*mfq_{j-1}^e_{j-1}*E_n -->...-->mfq_0^e_0*...*mfq_j^e_j*E_n
	if L_exp[0]>=0:
		L_chains_path.append(Chain_hor(osidh,0,L_chains_hor[0].j_center,L_chains_hor[0].L_plus[0:L_exp[0]],[]))
	else:
		L_chains_path.append(Chain_hor(osidh,0,L_chains_hor[0].j_center,[],L_chains_hor[0].L_minus[0:-L_exp[0]]))

	for j in range(1,t):
		chain=L_chains_hor[j]
		for k in range(j):
			chain=L_chains_path[k].action_chain(chain,L_exp[k],L_exp[j])
		L_chains_path.append(chain)
	return L_chains_path

def Action_hor(osidh,L_chains_hor,L_exp):
	r"""Given the horizontal chains at the bottom and the exponents, this function performs the
	action of the exponents on the central j-invariant of all the chains.

	INPUT:

	* osidh: OSIDH instanciation.

	* L_chains_hor: list of horizontal chains:
	mfq_j^-r*E_n-->...-->E_n -->...-->mfq_0^e_0*E_n

	* L_exp: list of exponents (integers in [-r,r])

	OUTPUT:

	Returns prod_{j=1}^t mfq_j^L_exp[j]*E_n.
	"""
	t=osidh.t
	L_chains_path=Action_hor_path(osidh,L_chains_hor,L_exp)
	if L_exp[t-1]==0:
		return L_chains_path[t-1].j_center
	if L_exp[t-1]>0:
		return L_chains_path[t-1].L_plus[L_exp[t-1]-1]
	else:
		return L_chains_path[t-1].L_minus[-L_exp[t-1]-1]


## Execution of the real OSIDH protocol
def OSIDH_exe(osidh,pub_chain):
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
	
	print("\nShared chains coincide: {0}".format(j_AB==j_BA))


## Protocol execution of the simplest version of OSIDH
def OSIDH_simple_exe(osidh,pub_chain):
	# Alice's secret key
	print("\nAlice:")
	L_exp_A=[randint(-osidh.r,osidh.r) for j in range(osidh.t)]
	print("Alice's secret key:")
	print(L_exp_A)
	chain_A=pub_chain.action(L_exp_A)	
	print("Alice's action on public chain complete.")

	# Bob's secret key
	print("\nBob:")
	L_exp_B=[randint(-osidh.r,osidh.r) for j in range(osidh.t)]
	print("Bob's secret key:")
	print(L_exp_B)
	chain_B=pub_chain.action(L_exp_B)	
	print("Bob's action on public chain complete.")

	# Alice's computation on Bob's chain
	chain_BA=chain_B.action(L_exp_A)
	print("\nAlice's action on Bob's chain complete.")

	# Bob's computation on Alice's chain
	chain_AB=chain_A.action(L_exp_B)
	print("\nBob's action on Alice's chain complete.")
	
	print("\nShared chains coincide: {0}".format(chain_AB.L_j==chain_BA.L_j))


## Obsolete functions: will be deleted soon
# Obsolete function: was used with bivariate polynomials
def add(L):
	r"""Intermediary function to add a list of polynomials and optimize addition cost 
	(getting very heavy when there are a lot of coefficients without optimization).

	INPUT: L: list of terms.

	OUTPUT: sum of the terms in L.
	"""

	k=len(L)
	if k==1:
		return L[0]
	else:
		return add(L[0:k//2])+add(L[k//2::])

def Phi(R,q,x,y):
	r"""
	Returns the modular polynomial of level q modulo p (implicit prime).

	INPUT:

	* R: bivariate polynomial ring above a finite field of characteristic p.

	* q: prime level

	* x, y: generators of R.

	OUTPUT: 

	Returns the modular polynomial of level q modulo p (implicit prime) in variables x and y, 
	retrieved from the database file phi_j_q.txt.
	"""

	try:
		with open(dirpath+"/Modular_polynomials/phi_j_{0}.txt".format(str(q)),"r",encoding="utf-8") as f:
			L_monom=monomials([x,y],[q+2,q+2])
			L_P=[]
			t_a=0
			t_b=0
			t_c=0
			t_d=0
			t_e=0
			for row in f:
				t1=time()
				L_row=row.split(" ")
				t2=time()
				L_ind=L_row[0][1:-1].split(",")
				t3=time()
				u,v=ZZ(L_ind[0]),ZZ(L_ind[1])
				t4=time()
				c=ZZ(L_row[1][0:-1])
				t5=time()
				if u!=v:
					L_P.append(c*(L_monom[(q+2)*u+v]+L_monom[(q+2)*v+u]))
				else:
					L_P.append(c*L_monom[(q+2)*u+v])
				t6=time()
				t_a+=t2-t1
				t_b+=t3-t2
				t_c+=t4-t3
				t_d+=t5-t4
				t_e+=t6-t5
			t1=time()
			P=add(L_P)
			t2=time()
			t_f=t2-t1
			return P,t_a,t_b,t_c,t_d,t_e,t_f
	except:
		raise ValueError("Prime q not in the database")

def Phi2(R,q,x,y):
	r"""
	Returns the modular polynomial of level q modulo p (implicit prime).

	INPUT:

	* R: bivariate polynomial ring above a finite field of characteristic p.

	* q: prime level

	* x, y: generators of R.

	OUTPUT: 

	Returns the modular polynomial of level q modulo p (implicit prime) in variables x and y, 
	retrieved from the database file phi_j_q.txt.
	"""

	try:
		with open(dirpath+"/Modular_polynomials/phi_{0}_modp.txt".format(str(q)),"r",encoding="utf-8") as f:
			P=R(0)
			i=0
			for row in f:
				if i%2==0:
					L_ind=row[1:-2].split(",")
					u,v=ZZ(L_ind[0]),ZZ(L_ind[1])
				else:
					c=ZZ(row[0:-1])
					if u!=v:
						P+=c*(x**u*y**v+x**v*y**u)
					else:
						P+=c*x**u*y**v
				i+=1
			return P
	except:
		raise ValueError("Prime q not in the database")

def RandomChain(R,S,j0,l,n):
	r"""
	Computes a descending l-isogeny chain of length n starting from j-invariant j0.

	INPUT:

	* R: bivariate polynomial ring above a finite field of characteristic p (base field).

	* S: univariate polynomial ring above the base field.

	* j0: element of the base field, starting j-invariant.

	* l: prime number.

	* n: positive integer.

	OUTPUT: 

	Returns the list of j-invariants representing the chain.
	"""

	L_j=[j0]

	x,y=R.gens()
	z=S.gen()
	phi_l=Phi(R,l,x,y)

	if n>=2:
		# Choice of the first j-invariant
		f=phi_l(j0,z)
		L_roots=f.roots()
		m=len(L_roots)
		i=randint(0,m-1)
		while L_roots[i][0]==j0:
			i=randint(0,m-1)
		L_j.append(L_roots[i][0])

		# Choice of the other j-invariant
		k=2
		while k<n:
			f=phi_l(L_j[k-1],z)
			L_roots=f.roots()
			m=len(L_roots)
			i=randint(0,m-1)
			while L_roots[i][0]==L_j[k-2]:
				i=randint(0,m-1)
			L_j.append(L_roots[i][0])
			k+=1
	return L_j



