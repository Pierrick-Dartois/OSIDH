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
	# Computation of sigma, the logarithm of the exponent in base p
	sigma=0
	for i in range(r):
		gi=B[i]**(p**sigma)
		while gi!=e:
			sigma+=1
			gi=gi**p
	print(sigma)
	C=[B]
	for i in range(sigma-1):
		C.append([g**p for g in C[-1]])
	print(C)
	H=[h]
	for j in range(1,sigma):
		H.append(H[j-1]**p)
	print(H)
	X=DL_exhaust(e,C[-1],H[-1],p,r)
	print(X)
	for j in range(1,sigma):
		hj=H[sigma-j-1]
		print(hj)
		print(C[sigma-j-1])
		for k in range(r):
			hj=hj*C[sigma-j-1][k]**(-X[k])
		print(hj)
		Y=DL_exhaust(e,C[-1],hj,p,r)
		print(Y)
		X=[X[i]+(p**j)*Y[i] for i in range(r)]
	return X

