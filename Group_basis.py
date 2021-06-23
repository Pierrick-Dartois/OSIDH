from sage.matrix.constructor import matrix
from sage.matrix.special import block_matrix
from sage.arith.all import lcm

### Discrete logarithm

## Echaustive search discrete logarithm in a p-group of exponent p
# Note that the result is in interval [0,p-1]**r
def DL_exhaust(e,B,h,p,r):
	if r==1:
		g0=B[0]
		g=e
		k=0
		while g!=h and k<p:
			g=g*g0
			k+=1
		if g==h:
			return [k%p]
		else:
			raise ValueError('Element g is not in the subgroup generated by the basis B')
	else:
		g=h
		C=B[0:r-1]
		gr=B[r-1]
		k=0
		while k<p:
			try:
				dl=DL_exhaust(e,C,g,p,r-1)
				return dl+[(p-k)%p] 
			except ValueError:
				g=gr*g
				k+=1
		raise ValueError('Element g is not in the subgroup generated by the basis B')


## Pohlig-Hellman in a p-group
def DL_PH(e,B,h,p,r):
	# Computation of the logarithm of the orders in base p and the basis of exponent p
	L_sigma=[0]*r
	C=B.copy()
	for i in range(r):
		gi=B[i]
		while gi!=e:
			L_sigma[i]+=1
			C[i]=gi
			gi=gi**p
	sigma=max(L_sigma)
	
	# Sorting the basis so that the orders decrease (using naive sorting algorithm)
	B1=B.copy()
	L_perm=list(range(r))# records the permutations
	for i in range(r):
		m=L_sigma[i]
		i_max=i
		for j in range(i,r):
			if L_sigma[j]>m:
				m=L_sigma[j]
				i_max=j
		L_sigma[i],L_sigma[i_max]=L_sigma[i_max],L_sigma[i]
		B1[i],B1[i_max]=B1[i_max],B1[i]
		C[i],C[i_max]=C[i_max],C[i]
		L_perm[i],L_perm[i_max]=L_perm[i_max],L_perm[i]	

	# Computing the powers of h
	H=[h]
	for i in range(1,sigma):
		H.append(H[i-1]**p)

	# Main Polling-Hellman loop
	X=[0]*r
	ri=0# Firts index such that L_sigma[ri]<=i
	for i in range(sigma-1,-1,-1):
		hi=H[i]
		while ri<r and L_sigma[ri]>=i+1:
			ri+=1
		for k in range(ri):
			hi=hi*B1[k]**(-(p**i)*X[k])
		Y=DL_exhaust(e,C[0:ri],hi,p,ri)
		for j in range(ri):
			X[j]+=p**(L_sigma[j]-i-1)*Y[j]
	
	# Inverse permutation 
	X_ret=[0]*r
	for i in range(r):
		X_ret[L_perm[i]]=X[i]
	return X_ret


### Basis of a group

## Basis of a p-group
def Basis_p(e,S,p):
	t=len(S)
	S1=S.copy()
	B=[]
	E_B=[]#log_p(orders) of B
	E=[0]*t#log_p(orders) of S
	for i in range(t):
		gi=S1[i]
		while gi!=e:
			E[i]+=1
			gi=gi**p
	E_max=max(E)

	# Adding first element in B
	if E_max==0:
		return B
	else:
		i0=E.index(E_max)
		B.append(S1[i0])
		E_B.append(E_max)
		del S1[i0]
		del E[i0]
		t=t-1

	# Main loop with orthogonalization process
	while t>0:
		r=len(B)
		S_update=[]
		E_update=[]
		for i in range(t):
			if E[i]>0:
				b_continue=True
				e_i=0
				C=B.copy()
				hi=S1[i]
				while b_continue and e_i<E[i]:
					ri=0
					while ri<r and E_B[ri]>=e_i+1:
						ri+=1
					try:
						Xi=DL_PH(e,C[0:ri],hi,p,ri)
						b_continue=False
					except:
						C=[c**p for c in C]
						hi=hi**p
						e_i+=1
				s=S1[i]
				if e_i<E[i]:
					for j in range(ri):
						s*=B[j]**(-Xi[j])
				S_update.append(s)
				E_update.append(e_i)
		E=E_update.copy()
		S1=S_update.copy()
		# Updating S and B as previously
		if len(E)==0:
			return B
		else:
			E_max=max(E)
			if E_max==0:
				return B
			else:
				i0=E.index(E_max)
				B.append(S1[i0])
				E_B.append(E_max)
				del S1[i0]
				del E[i0]
				t=len(S1)
	return B


## Testing function for Basis_p
def IsBasis(e,B,C,p):
	r=len(B)
	if len(C)!=r:
		return False

	# Computing the log in base p of the orders of B and C
	L_sigma_B=[0]*r
	L_sigma_C=[0]*r
	for i in range(r):
		gi=B[i]
		while gi!=e:
			L_sigma_B[i]+=1
			gi=gi**p
	for i in range(r):
		gi=C[i]
		while gi!=e:
			L_sigma_C[i]+=1
			gi=gi**p
	# We should have L_sigma_B.sort()=L_sigma_C.sort() by unicity of the invariant factors
	if L_sigma_B.sort()!=L_sigma_C.sort():
		return False

	# Tests if all elements of B are generated by C and conclude by cardinality
	for i in range(r):
		try:
			X=DL_PH(e,C,B[i],p,r)
		except ValueError:
			return False
	return True


## 



### Get relations from a basis

## Returns the basis of the lattice of vectors e such that M_i e \equiv B_i for all i
def Lattice_basis(M,b):
	r=len(b)
	t=M.ncols()
	L_A=[]
	for i in range(r):
		Ai=block_matrix([[matrix(M[i])],[b[i]]])
		# Here, the generator vectors are the lines and not the columns (because of HNF)
		L_A.append(Ai)
	m=lcm(b)
	A=block_matrix([[(m//b[i])*L_A[i]] for i in range(r)])
	A1=A.hermite_form()
	return m*A1[0:t,0:t].inverse()


## Testing function for Lattice_basis
def IsGoodLatticeBasis(B,M,b):
	r=len(b)
	t=M.ncols()
	m=1
	for i in range(r):
		m*=b[i]
	if m!=B.det():
		return False
	C=M*B
	for i in range(r):
		for j in range(t):
			if C[i,j]%b[i]!=0:
				return False
	return True






