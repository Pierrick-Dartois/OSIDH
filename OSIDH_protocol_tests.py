import sys
sys.path.append("Documents/Codes/OSIDH")
from sage.all import *
from OSIDH_protocol import *
from time import time
import pickle

dirpath="Documents/Codes/OSIDH"

### Usual toy parameters 
n=28
t=10
l=2
r=3
# so that h\simeq l**n\simeq (2r+1)**t
d_K=-4
print("Testing with usual toy parameters: n={0}, t={1}, l={2}, r={3}, d_K={4}".format(n,t,l,r,d_K))

## Class OSIDH
# Test 1: test __init__
try:
	t1=time()
	osidh=OSIDH(n,t,l,r,d_K)
	t2=time()
	print("Test 1 (creation OSIDH): success")
	print("Time test 1: {0} s".format(t2-t1))
except:
	print("Test 1 (creation OSIDH): failure")

# Test 1-bis: pickling
try:
	t1=time()
	osidh.save(dirpath+"/pickle_OSIDH.txt")
	t2=time()
	print("Test 1-bis (pickling): success")
	print("Time test 1-bis: {0} s".format(t2-t1))
except:
	print("Test 1-bis (pickling): failure")

# Test 1-ter: unpickling
try:
	t1=time()
	with open(dirpath+"/pickle_OSIDH.txt","rb") as f:
		osidh=pickle.load(f)
	t2=time()
	print("Test 1-ter (unpickling): success")
	print("Time test 1-ter: {0} s".format(t2-t1))
except:
	print("Test 1-ter (unpickling): failure")

## Class Chain
# Test 2: chain creation
try:
	t1=time()
	C=Chain(osidh)
	t2=time()
	print("Test 2 (chain creation): success")
	print("Time test 2: {0} s".format(t2-t1))
except:
	print("Test 2 (chain creation): failure")

## q-action
# Test 3: action on the chian by the smallest prime
try:
	t1=time()
	D=C.action_prime(0)
	t2=time()
	print("Test 3 (q-action): success")
	print("Time test 3: {0} s".format(t2-t1))
except:
	print("Test 3 (q-action): failure")

# Test 4: action on the chian by the biggest prime
try:
	t1=time()
	D=C.action_prime(t-1)
	t2=time()
	print("Test 4 (q-action): success")
	print("Time test 4: {0} s".format(t2-t1))
except:
	print("Test 4 (q-action): failure")


### Toy parameters with different l
n=18
t=10
l=3
r=3
d_K=-4
print("\nTesting with different l: n={0}, t={1}, l={2}, r={3}, d_K={4}".format(n,t,l,r,d_K))

## Class OSIDH
# Test 1: test __init__
try:
	t1=time()
	osidh=OSIDH(n,t,l,r,d_K)
	t2=time()
	print("Test 1 (creation OSIDH): success")
	print("Time test 1: {0} s".format(t2-t1))
except:
	print("Test 1 (creation OSIDH): failure")

## Class Chain
# Test 2: chain creation
try:
	t1=time()
	C=Chain(osidh)
	t2=time()
	print("Test 2 (chain creation): success")
	print("Time test 2: {0} s".format(t2-t1))
except:
	print("Test 2 (chain creation): failure")

## q-action
# Test 3: action on the chian by the smallest prime
try:
	t1=time()
	D=C.action_prime(0)
	t2=time()
	print("Test 3 (q-action): success")
	print("Time test 3: {0} s".format(t2-t1))
except:
	print("Test 3 (q-action): failure")

### Toy parameters with different d_K
n=28
t=10
l=2
r=3
d_K=-3
print("\nTesting with different d_K: n={0}, t={1}, l={2}, r={3}, d_K={4}".format(n,t,l,r,d_K))

## Class OSIDH
# Test 1: test __init__
try:
	t1=time()
	osidh=OSIDH(n,t,l,r,d_K)
	t2=time()
	print("Test 1 (creation OSIDH): success")
	print("Time test 1: {0} s".format(t2-t1))
except:
	print("Test 1 (creation OSIDH): failure")

## Class Chain
# Test 2: chain creation
try:
	t1=time()
	C=Chain(osidh)
	t2=time()
	print("Test 2 (chain creation): success")
	print("Time test 2: {0} s".format(t2-t1))
except:
	print("Test 2 (chain creation): failure")

## q-action
# Test 3: action on the chian by the smallest prime
try:
	t1=time()
	D=C.action_prime(0)
	t2=time()
	print("Test 3 (q-action): success")
	print("Time test 3: {0} s".format(t2-t1))
except:
	print("Test 3 (q-action): failure")

### Toy parameters with different d_K and l
n=18
t=10
l=3
r=3
d_K=-3
print("\nTesting with different d_K and l: n={0}, t={1}, l={2}, r={3}, d_K={4}".format(n,t,l,r,d_K))

## Class OSIDH
# Test 1: test __init__
try:
	t1=time()
	osidh=OSIDH(n,t,l,r,d_K)
	t2=time()
	print("Test 1 (creation OSIDH): success")
	print("Time test 1: {0} s".format(t2-t1))
except:
	print("Test 1 (creation OSIDH): failure")

## Class Chain
# Test 2: chain creation
try:
	t1=time()
	C=Chain(osidh)
	t2=time()
	print("Test 2 (chain creation): success")
	print("Time test 2: {0} s".format(t2-t1))
except:
	print("Test 2 (chain creation): failure")

## q-action
# Test 3: action on the chian by the smallest prime
try:
	t1=time()
	D=C.action_prime(0)
	t2=time()
	print("Test 3 (q-action): success")
	print("Time test 3: {0} s".format(t2-t1))
except:
	print("Test 3 (q-action): failure")

### Testing the whole protocol
### Usual toy parameters 
n=28
t=10
l=2
r=3
d_K=-4
print("\nTesting OSIDH protocol with: n={0}, t={1}, l={2}, r={3}, d_K={4}".format(n,t,l,r,d_K))
osidh=OSIDH(n,t,l,r,d_K)
pub_chain=Chain(osidh)

t1=time()
OSIDH_exe(osidh,pub_chain)
t2=time()

print("OSIDH execution time: {0} s".format(t2-t1))
