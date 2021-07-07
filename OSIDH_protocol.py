from sage.all import *

dirpath="Documents/Codes/OSIDH"

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

def Ladder_q(R,S,L_chain,l,q):
	pass



