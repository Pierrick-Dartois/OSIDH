import sys
sys.path.append("Documents/Codes/OSIDH")
from sage.all import *
from OSIDH_protocol import *
from time import time
import pickle

dirpath="Documents/Codes/OSIDH"


## Test __init__
print("Testing constructor")

# Test 1:
# Usual toy Parameters
n=28
t=10
l=2
r=3
# so that h\simeq l**n\simeq (2r+1)**t
d_K=-4

try:
	t1=time()
	#prot=OSIDH(n,t,l,r,d_K)
	t2=time()
	print("Test 1: success")
	print("Time test 1: {0} s".format(t2-t1))
except:
	print("Test 1: failure")

# Test 1-bis: pickling
try:
	t1=time()
	#prot.save(dirpath+"/pickle_OSIDH.txt")
	t2=time()
	print("Test 1-bis (pickling): success")
	print("Time test 1-bis: {0} s".format(t2-t1))
except:
	print("Test 1-bis (pickling): failure")

# Test 1-ter: unpickling
try:
	t1=time()
	#with open(dirpath+"/pickle_OSIDH.txt","rb") as f:
		#prot2=pickle.load(f)
	t2=time()
	print("Test 1-ter (unpickling): success")
	print("Time test 1-ter: {0} s".format(t2-t1))
except:
	print("Test 1-ter (unpickling): failure")


# Test 2:
# Toy Parameters with different l
n=18
t=10
l=3
r=3
# so that h\simeq l**n\simeq (2r+1)**t
d_K=-4

try:
	t1=time()
	#prot=OSIDH(n,t,l,r,d_K)
	t2=time()
	print("Test 2: success")
	print("Time test 2: {0} s".format(t2-t1))
except:
	print("Test 2: failure")


# Test 3:
# Toy Parameters with different d_K
n=28
t=10
l=2
r=3
# so that h\simeq l**n\simeq (2r+1)**t
d_K=-3

try:
	t1=time()
	#prot=OSIDH(n,t,l,r,d_K)
	t2=time()
	print("Test 3: success")
	print("Time test 3: {0} s".format(t2-t1))
except:
	print("Test 3: failure")

# Test 4:
# Toy Parameters with different d_K and l
n=18
t=10
l=3
r=3
# so that h\simeq l**n\simeq (2r+1)**t
d_K=-3

try:
	t1=time()
	#prot=OSIDH(n,t,l,r,d_K)
	t2=time()
	print("Test 4: success")
	print("Time test 4: {0} s".format(t2-t1))
except:
	print("Test 4: failure")

# Test 5:
# Big polynomial
n=256
t=74
l=2
r=5
d_K=-4

L_q=[]
q=2
prod=l
for j in range(t):
	while kronecker(d_K,q)==-1 or q==l:
		q=next_prime(q)
	L_q.append(q)
	prod*=q
	q=next_prime(q)

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

F=GF(p**2,"a",proof="False")
Fxy=PolynomialRing(F,["x","y"])
x,y=Fxy.gens()

t1=time()
P,t_a,t_b,t_c,t_d,t_e,t_f=Phi(Fxy,881,x,y)
t2=time()

print("Test 5: {0} s".format(t2-t1))

