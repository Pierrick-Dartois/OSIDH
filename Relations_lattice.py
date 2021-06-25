import sys
sys.path.append('Documents/Codes/OSIDH')
from Group_basis import *
from Class_group import *
from sage.all import *

n=100
t=50
l=2
d_K=-4
f=l**n
d=d_K*f**2

L_q=[]
q=3
for j in range(t):
	while not Mod(d_K,q).is_square() or q==l:
		q=next_prime(q)
	L_q.append(q)
	q=next_prime(q)

Q=[IdealClass(d,q) for q in L_q]
for mfq in Q:
	print(str(mfq))

