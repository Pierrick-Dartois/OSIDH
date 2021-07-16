import sys
sys.path.append("Documents/Codes/OSIDH")
from sage.all import *
from OSIDH_protocol import *
from time import time

# Usual toy Parameters
n=28
t=10
l=2
r=3
# so that h\simeq l**n\simeq (2r+1)**t
d_K=-4

osidh=OSIDH(n,t,l,r,d_K)
C=Chain(osidh)
mfq=gp.qfbprimeform(d_K,osidh.L_q[t-1])
t1=time()
j=C.action_torsion(mfq,t-1,5)
t2=time()
print(t2-t1)