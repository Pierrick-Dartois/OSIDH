from sage.all import *

dirpath="Documents/Codes/OSIDH"

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
		
		if not l.is_prime():
			raise ValueError("l should be a prime")
		if d_K not in [3,4]:
			raise ValueError("d_K should be -3 or -4")

		self.n=n
		self.t=t
		self.l=l
		self.r=r
		self.d_K=d_K

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

		# Field of definition and polynomial rings
		self.F=GF(p**2,"a",proof="False")
		self.Fxy=PolynomialRing(F,["x","y"])
		self.Fz=PolynomialRing(F,"z")

		# Modular polynomials
		










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
		with open(dirpath+"/phi_j_upto_300/phi_j_{0}.txt".format(str(q)),"r",encoding="utf-8") as f:
			P=R(0)
			for row in f:
				L_row=row.split(" ")
				L_ind=L_row[0][1:-1].split(",")
				u,v=ZZ(L_ind[0]),ZZ(L_ind[1])
				c=R(L_row[1][0:-1])
				if u!=v:
					P+=c*(x**u*y**v+x**v*y**u)
				else:
					P+=c*x**u*y**v
			return P
	except:
		raise ValueError("Prime q not in the database")

def Chain(R,S,j0,l,n):
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

def IsogenyFromCurves(E0,E1,l):
	pass


def Ladder_q(R,S,L_chain,l,q):
	pass



