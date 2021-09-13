from sage.all import *
from time import time
from OSIDH_protocol import *
from OSIDH_attack import *

## Non cyclic case
# Parameters
n=10
t=10
l=2
r=3
d_K=-3
print("\nTest non-cyclic case: l={0}, d_K={1}\n".format(l,d_K))

osidh=OSIDH(n,t,l,r,d_K)
for i in range(10):
	pub_chain=Chain(osidh)
	L_exp=[randint(-r,r) for j in range(t)]
	chain_ex=pub_chain.action(L_exp)
	L_exp_attack=attack_chain(pub_chain,chain_ex)
