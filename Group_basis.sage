### Discrete logarithm
## Echaustive search discrete logarithm in a p-group of exponent p
# Note that the result is in interval [1,p]^r, instead of conventionnal [0,p-1]^r
def DL_exhaust(B,h,p,r):
	if r==1:
		g0=B[0]
		g=g0
		k=1
		while g!=h and k<p+1:
			g=g*g0
			k+=1
		if g==h:
			return [k]
		else:
			raise ValueError('Element g is not in the subgroup generated by the basis B')
	else:
		g=h
		C=B[0:r-1]
		gr=B[r-1]
		k=0
		while k<p:
			try:
				dl=DL_exhaust(C,g,p,r-1)
				return dl+[p-k] 
			except ValueError:
				g=gr*g
				k+=1
		raise ValueError('Element g is not in the subgroup generated by the basis B')

###Tests
##Discrete logarithm
#DL_exhaust

#Test 1: A simple cyclic group
G=CyclicPermutationGroup(23)
g=G.gen()
h=G.random_element()
k=DL_exhaust([g],h,23,1)
if g^(k[0])==h:
	print("Test 1: success")
else:
	print("Test 1: failure")

#Test 2: for r=2
G1=CyclicPermutationGroup(19)
G=cartesian_product([G1,G1])
g0=G1.gen()
B=[G((g0,G1.one())),G((G1.one(),g0))]
h=G.random_element()
k=DL_exhaust(B,h,19,2)
if B[0]^(k[0])*B[1]^(k[1])==h:
	print("Test 2: success")
else:
	print("Test 2: failure")

#Test 3: for r=3
G1=CyclicPermutationGroup(31)
G=cartesian_product([G1,G1,G1])
g0=G1.gen()
e=G1.one()
B=[G((g0,e,e)),G((e,g0,e)),G((e,e,g0))]
h=G.random_element()
k=DL_exhaust(B,h,31,3)
if B[0]^(k[0])*B[1]^(k[1])*B[2]^(k[2])==h:
	print("Test 3: success")
else:
	print("Test 3: failure")

#Test 4: in bigger groups


