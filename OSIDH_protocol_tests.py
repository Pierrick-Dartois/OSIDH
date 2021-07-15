import sys
sys.path.append("Documents/Codes/OSIDH")
from sage.all import *
from OSIDH_protocol import *
from time import time
import pickle

dirpath="Documents/Codes/OSIDH"

### Class OSIDH

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
	prot=OSIDH(n,t,l,r,d_K)
	t2=time()
	print("Test 1: success")
	print("Time test 1: {0} s".format(t2-t1))
except:
	print("Test 1: failure")

# Test 1-bis: pickling
try:
	t1=time()
	prot.save(dirpath+"/pickle_OSIDH.txt")
	t2=time()
	print("Test 1-bis (pickling): success")
	print("Time test 1-bis: {0} s".format(t2-t1))
except:
	print("Test 1-bis (pickling): failure")

# Test 1-ter: unpickling
try:
	t1=time()
	with open(dirpath+"/pickle_OSIDH.txt","rb") as f:
		prot2=pickle.load(f)
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
	prot=OSIDH(n,t,l,r,d_K)
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
	prot=OSIDH(n,t,l,r,d_K)
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
	prot=OSIDH(n,t,l,r,d_K)
	t2=time()
	print("Test 4: success")
	print("Time test 4: {0} s".format(t2-t1))
except:
	print("Test 4: failure")


### Class Chain

print("\nTesting class Chain")

## Test 1: chain creation
# Usual toy Parameters
n=28
t=10
l=2
r=3
d_K=-4

osidh=OSIDH(n,t,l,r,d_K)

try:
	t1=time()
	C=Chain(osidh)
	t2=time()
	print("Test 1: success")
	print("Time test 1: {0} s".format(t2-t1))
except:
	print("Test 1: failure")

mfq=gp.qfbprimeform(d_K,osidh.L_q[t-1])
t1=time()
#D=C.action_prime(mfq,t-1)
#j=C.action_torsion(mfq,1)
t2=time()
print(t2-t1)