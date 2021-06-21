import sys
sys.path.append('Documents/Codes/OSIDH')
from Group_basis import *

##Discrete logarithm
#DL_exhaust
print("Testing DL_exhaust")
#Test 1: A simple cyclic group
G=CyclicPermutationGroup(23)
g=G.gen()
h=G.random_element()
k=DL_exhaust(G.one(),[g],h,23,1)
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
k=DL_exhaust(G.one(),B,h,19,2)
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
k=DL_exhaust(G.one(),B,h,31,3)
if B[0]^(k[0])*B[1]^(k[1])*B[2]^(k[2])==h:
	print("Test 3: success")
else:
	print("Test 3: failure")

#Test 4: in bigger groups
G1=CyclicPermutationGroup(5)
G2=CyclicPermutationGroup(5^2)
g1=G1.gen()
e1=G1.one()
g2=G2.gen()
e2=G2.one()
G=cartesian_product([G1,G2])
B=[G((g1,e2)),G((e1,g2^5))]#g2^5 keeps exponent 5
h=G((G1.random_element(),G2.random_element()^5))
k=DL_exhaust(G.one(),B,h,5,2)
if B[0]^(k[0])*B[1]^(k[1])==h:
	print("Test 4: success")
else:
	print("Test 4: failure")

#Test 5: raise an exception for r=1
G=CyclicPermutationGroup(23^2)
g=G.gen()
h=g^21
try :
	DL_exhaust(G.one(),[g^23],h,23,1)
	print("Test 5: failure")
except ValueError:
	print("Test 5: success")

#Test 6: raise an exception for r=2
G1=CyclicPermutationGroup(23^2)
g0=G1.gen()
e=G1.one()
G=cartesian_product([G1,G1])
B=[G((g0^23,e)),G((e,g0^23))]
h=G((g0^17,g0^25))
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
if B[0]^(k[0])*B[1]^(k[1])==h:
	print("Test 7: success")
else:
	print("Test 7: failure")


##DL_PH
print("Testing DL_PH")

#Test 1: r=1
G=AbelianGroup([3^2])
g=G.gen()
h=G.random_element()
k=DL_PH(G.one(),[g],h,3,1)
if g^(k[0])==h:
	print("Test 1: success")
else:
	print("Test 1: failure")

#Test 2: r=2
G=AbelianGroup([3^2,3^3])
(g1,g2)=G.gens()
B=[g1,g2]
e=G.one()
h=G.random_element()
k=DL_PH(e,B,h,3,2)
if B[0]^(k[0])*B[1]^(k[1])==h:
	print("Test 2: success")
else:
	print("Test 2: failure")

#Test 3: r=5, big exponents
G=AbelianGroup([3^7,3^2,3^5,3^4,3^8])
B=list(G.gens())
e=G.one()
h=G.random_element()
k=DL_PH(e,B,h,3,5)
res=e
for j in range(5):
	res*=B[j]^k[j]
if res==h:
	print("Test 3: success")
else:
	print("Test 3: failure")

#Test 4: p=2, r=5, big exponents
G=AbelianGroup([2^17,2^12,2^15,2^24,2^18])
B=list(G.gens())
e=G.one()
h=G.random_element()
k=DL_PH(e,B,h,2,5)
res=e
for j in range(5):
	res*=B[j]^k[j]
if res==h:
	print("Test 4: success")
else:
	print("Test 4: failure")
