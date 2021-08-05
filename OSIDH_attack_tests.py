import sys
sys.path.append("Documents/Codes/OSIDH")
from sage.all import *
from time import time
from OSIDH_protocol import *
from OSIDH_attack import *

# Parameters
n=28
t=10
l=2
r=3
d_K=-4
d=d_K*l**(2*n)

osidh=OSIDH(n,t,l,r,d_K)
pub_chain=Chain(osidh)
OSIDH_attack_exe(osidh,pub_chain)