import sys
sys.path.append("Documents/Codes/OSIDH")
from sage.all import *
from time import time
from OSIDH_protocol import *
from OSIDH_attack import *

## Cyclic case
# Parameters
n=15# To shorten the attack
t=10
l=2
r=3
d_K=-4
print("Test cyclic case: l={0}, d_K={1}".format(l,d_K))

osidh=OSIDH(n,t,l,r,d_K)
pub_chain=Chain(osidh)
L_exp=[randint(-osidh.r,osidh.r) for j in range(osidh.t)]
chain_ex=pub_chain.action(L_exp)
L_exp_attack=attack_chain(pub_chain,chain_ex)
chain_ex_attack=pub_chain.action(L_exp_attack)
print("Attack is correct: {0}".format(chain_ex_attack.L_j==chain_ex.L_j))

## Non cyclic case
# Parameters
n=15# To shorten the attack
t=10
l=2
r=3
d_K=-3
print("Test non-cyclic case: l={0}, d_K={1}".format(l,d_K))

osidh=OSIDH(n,t,l,r,d_K)
pub_chain=Chain(osidh)
L_exp=[randint(-osidh.r,osidh.r) for j in range(osidh.t)]
chain_ex=pub_chain.action(L_exp)
L_exp_attack=attack_chain(pub_chain,chain_ex)
chain_ex_attack=pub_chain.action(L_exp_attack)
print("Attack is correct: {0}".format(chain_ex_attack.L_j==chain_ex.L_j))

## Cyclic case, different l
# Parameters
n=15# To shorten the attack
t=10
l=3
r=3
d_K=-3
print("Cyclic case: l={0}, d_K={1}".format(l,d_K))

osidh=OSIDH(n,t,l,r,d_K)
pub_chain=Chain(osidh)
L_exp=[randint(-osidh.r,osidh.r) for j in range(osidh.t)]
chain_ex=pub_chain.action(L_exp)
L_exp_attack=attack_chain(pub_chain,chain_ex)
chain_ex_attack=pub_chain.action(L_exp_attack)
print("Attack is correct: {0}".format(chain_ex_attack.L_j==chain_ex.L_j))
