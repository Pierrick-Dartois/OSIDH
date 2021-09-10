import sys
sys.path.append('Documents/Codes/OSIDH')
from Group_basis import *
from sage.all import *

##Discrete logarithm
#DL_exhaust
print("Testing DL_exhaust (test 7 may raise a ValueError)")

#Test 1: A simple cyclic group
G=CyclicPermutationGroup(23)
g=G.gen()
h=G.random_element()
k=DL_exhaust(G.one(),[g],h,23,1)
if g**(k[0])==h:
	print("Test 1: success")
else:
	print("Test 1: failure")

#Test 2: for r=2
G1=CyclicPermutationGroup(19)
G=cartesian_product([G1,G1])
g0=G1.gen()
B=[G((g0,G1.one())),G((G1.one(),g0))]
h=G.random_element()
k=DL_exhaust(G.one(),B,h,19,2)
if B[0]**(k[0])*B[1]**(k[1])==h:
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
k=DL_exhaust(G.one(),B,h,31,3)
if B[0]**(k[0])*B[1]**(k[1])*B[2]**(k[2])==h:
	print("Test 3: success")
else:
	print("Test 3: failure")

#Test 4: in bigger groups
G1=CyclicPermutationGroup(5)
G2=CyclicPermutationGroup(5**2)
g1=G1.gen()
e1=G1.one()
g2=G2.gen()
e2=G2.one()
G=cartesian_product([G1,G2])
B=[G((g1,e2)),G((e1,g2**5))]#g2**5 keeps exponent 5
h=G((G1.random_element(),G2.random_element()**5))
k=DL_exhaust(G.one(),B,h,5,2)
if B[0]**(k[0])*B[1]**(k[1])==h:
	print("Test 4: success")
else:
	print("Test 4: failure")

#Test 5: raise an exception for r=1
G=CyclicPermutationGroup(23**2)
g=G.gen()
h=g**21
try :
	DL_exhaust(G.one(),[g**23],h,23,1)
	print("Test 5: failure")
except ValueError:
	print("Test 5: success")

#Test 6: raise an exception for r=2
G1=CyclicPermutationGroup(23**2)
g0=G1.gen()
e=G1.one()
G=cartesian_product([G1,G1])
B=[G((g0**23,e)),G((e,g0**23))]
h=G((g0**17,g0**25))
try :
	DL_exhaust(G.one(),B,h,23,2)
	print("Test 6: failure")
except ValueError:
	print("Test 6: success")

#Test 7: mixing elements (there is a low probability of failure)
G1=CyclicPermutationGroup(17)
g0=G1.gen()
e=G1.one()
G=cartesian_product([G1,G1])
B=[G.random_element(),G.random_element()]
h=G.random_element()
k=DL_exhaust(G.one(),B,h,17,2)
if B[0]**(k[0])*B[1]**(k[1])==h:
	print("Test 7: success")
else:
	print("Test 7: failure")


##DL_PH
print("Testing DL_PH")

#Test 1: r=1
G=AbelianGroup([3**2])
g=G.gen()
h=G.random_element()
k=DL_PH(G.one(),[g],h,3,1)
if g**(k[0])==h:
	print("Test 1: success")
else:
	print("Test 1: failure")

#Test 2: r=2
G=AbelianGroup([3**2,3**3])
(g1,g2)=G.gens()
B=[g1,g2]
e=G.one()
h=G.random_element()
k=DL_PH(e,B,h,3,2)
if B[0]**(k[0])*B[1]**(k[1])==h:
	print("Test 2: success")
else:
	print("Test 2: failure")

#Test 3: r=5, big exponents
G=AbelianGroup([3**7,3**2,3**5,3**4,3**8])
B=list(G.gens())
e=G.one()
h=G.random_element()
k=DL_PH(e,B,h,3,5)
res=e
for j in range(5):
	res*=B[j]**k[j]
if res==h:
	print("Test 3: success")
else:
	print("Test 3: failure")

#Test 4: p=2, r=5, big exponents
G=AbelianGroup([2**17,2**12,2**15,2**24,2**18])
B=list(G.gens())
e=G.one()
h=G.random_element()
k=DL_PH(e,B,h,2,5)
res=e
for j in range(5):
	res*=B[j]**k[j]
if res==h:
	print("Test 4: success")
else:
	print("Test 4: failure")

##DL
print("Testing DL")

# Test 1: r=1, 2 factors
G=AbelianGroup([2**5*3**4])
g=G.gen()
h=G.random_element()
k=DL(G.one(),[g],h,[(2,5),(3,4)],1)
if g**(k[0])==h:
	print("Test 1: success")
else:
	print("Test 1: failure")

# Test 2: r=2, 2 factors
G=AbelianGroup([2**5*3**4,2**2*3**3])
g0,g1=G.gens()
h=G.random_element()
k=DL(G.one(),[g0,g1],h,[(2,7),(3,7)],2)
if g0**(k[0])*g1**(k[1])==h:
	print("Test 2: success")
else:
	print("Test 2: failure")

# Test 3: r=2, 4 factors
G=AbelianGroup([2**5*3**4*7**5,2**2*3**3*13**14])
g0,g1=G.gens()
h=G.random_element()
k=DL(G.one(),[g0,g1],h,[(2,7),(3,7),(7,5),(13,14)],2)
if g0**(k[0])*g1**(k[1])==h:
	print("Test 3: success")
else:
	print("Test 3: failure")

# Test 4: r=1, 4 factors
G=AbelianGroup([2**7*3**7*7**5*13**14])
g=G.gen()
h=G.random_element()
k=DL(G.one(),[g],h,[(2,7),(3,7),(7,5),(13,14)],1)
if g**(k[0])==h:
	print("Test 4: success")
else:
	print("Test 4: failure")


## Basis_p
print("Testing Basis_p: the tests are probabilistic so they may fail."
	+"\n See test 6 for multiple testing taking into account"+
	"\n this probability of failure (should succeed).")

# Test 1 : r=2, t=3
G=AbelianGroup([3**2,3**3])
e=G.one()
S=[G.random_element() for i in range(3)]
B=list(G.gens())
C,E_C=Basis_p(e,S,3)
if IsBasis(e,B,C,E_C,3):
	print("Test 1: success")
else:
	print("Test 1: failure")

# Test 2: r=5, t=10
G=AbelianGroup([3**7,3**2,3**5,3**4,3**8])
e=G.one()
S=[G.random_element() for i in range(10)]
B=list(G.gens())
C,E_C=Basis_p(e,S,3)
if IsBasis(e,B,C,E_C,3):
	print("Test 2: success")
else:
	print("Test 2: failure")

# Test 3: r=1, t=1
G=AbelianGroup([3**7])
e=G.one()
S=[G.random_element()]
B=list(G.gens())
C,E_C=Basis_p(e,S,3)
if IsBasis(e,B,C,E_C,3):
	print("Test 3: success")
else:
	print("Test 3: failure")

# Test 4: r=1, t=2
G=AbelianGroup([3**7])
e=G.one()
S=[G.random_element(),G.random_element()]
B=list(G.gens())
C,E_C=Basis_p(e,S,3)
if IsBasis(e,B,C,E_C,3):
	print("Test 4: success")
else:
	print("Test 4: failure")

# Test 5 : Realistic test for ideal class groups : p=2, r=2 big cyclic factor * Z/2Z
G=AbelianGroup([2**100,2])
e=G.one()
S=[G.random_element() for i in range(50)]
B=list(G.gens())
C,E_C=Basis_p(e,S,2)
if IsBasis(e,B,C,E_C,2):
	print("Test 5: success")
else:
	print("Test 5: failure")

# Test 6 : Statistical test, to see if we get a basis when the initial set is a basis
p=5
n=10
G=AbelianGroup([p**n,p**n])
e=G.one()
B=list(G.gens())
b_test=True
for i in range(100):
	S=[G.random_element() for i in range(2)]
	C,E_C=Basis_p(e,S,p)
	if IsBasis(e,B,C,E_C,p)!=IsBasis2(e,B,S,p):
		b_test=False
if b_test:
	print("Test 6: success")
else:
	print("Test 6: failure")


## Basis
print("Testing Basis")

# Test 1:
G=AbelianGroup([2**20,2**5,3**25,3**50,5**6])
e=G.one()
S=[G.random_element() for i in range(50)]
N=(2**20)*(2**5)*(3**25)*(3**50)*(5**6)
L_factors=[(2,25),(3,75),(5,6)]
B,L_orders=Basis(e,S,N,L_factors)
if L_orders==[(2,20),(2,5),(3,50),(3,25),(5,6)]:
	print("Test 1: success")
else:
	print("Test 1: failure")


## DL_matrix_sylow
print("Testing DL_matrix_sylow")

# Test 1: uses Test 1 of Basis
M=DL_matrix_sylow(e,S,B,L_orders)
if IsGoodDL_matrix(e,B,M,S):
	print("Test 1: success")
else:
	print("Test 1: failure")

## DL_matrix
print("Testing DL_matrix and Basis_clg")

# Test 1: uses Test 1 of Basis
M=DL_matrix(e,S,B,L_factors)
if IsGoodDL_matrix(e,B,M,S):
	print("Test 1 (DL_matrix): success")
else:
	print("Test 1 (DL_matrix): failure")

# Test 2: with Basis_clg Cl(O_28) for d_K=-11
G=AbelianGroup([2,3*2**26])
L_factors=[(3,1),(2,27)]
e=G.one()
S=[G.random_element() for i in range(50)]
# Find a generator of the second component Z/(3*2**26)Z
B=Basis_clg(e,S,2,28,3)
M=DL_matrix(e,S,B,L_factors)
if IsGoodDL_matrix(e,B,M,S):
	print("Test 2 (DL_matrix and Basis_clg): success")
else:
	print("Test 2 (DL_matrix and Basis_clg): failure")

## Lattice basis
print("Testing lattice basis (test 2 may fail)")

# Test 1 : on the group of test 2 of DL_matrix
b=[3*2**26,2]
B=Lattice_basis(M,b)
if IsGoodLatticeBasis(B,M,b):
	print("Test 1: success")
else:
	print("Test 1: failure")

# Test 2 : multiple examples 
b_test=True
for i in range(100):
	r=randint(1,5)
	t=randint(2*r,30)
	M=random_matrix(ZZ,r,t)
	p=random_prime(100)
	q=random_prime(100)
	b=[p**(randint(0,r))*q**(randint(0,r)) for i in range(r)]
	B=Lattice_basis(M,b)
	if not IsGoodLatticeBasis(B,M,b):
		b_test=False
if b_test:
	print("Test 2: success")
else:
	print("Test 2: failure")



