from sage.all import *
from time import time
from OSIDH_protocol import *
from OSIDH_attack import *

## Cyclic case
# Parameters
n=28
t=10
l=2
r=3
d_K=-4
print("Test cyclic case: l={0}, d_K={1}\n".format(l,d_K))

osidh=OSIDH(n,t,l,r,d_K)
pub_chain=Chain(osidh)
OSIDH_attack_exe(osidh,pub_chain)

## Non cyclic case
# Parameters
n=28
t=10
l=2
r=3
d_K=-3
print("\nTest non-cyclic case: l={0}, d_K={1}\n".format(l,d_K))

osidh=OSIDH(n,t,l,r,d_K)
pub_chain=Chain(osidh)
OSIDH_attack_exe(osidh,pub_chain)

## Cyclic case, different l (d_K=-3)
# Parameters
n=18
t=10
l=3
r=3
d_K=-3
print("\nCyclic case: l={0}, d_K={1}\n".format(l,d_K))

osidh=OSIDH(n,t,l,r,d_K)
pub_chain=Chain(osidh)
OSIDH_attack_exe(osidh,pub_chain)

## Cyclic case, different l (d_K=-4)
# Parameters
n=12
t=10
l=5
r=3
d_K=-4
print("\nCyclic case: l={0}, d_K={1}\n".format(l,d_K))

osidh=OSIDH(n,t,l,r,d_K)
pub_chain=Chain(osidh)
OSIDH_attack_exe(osidh,pub_chain)
