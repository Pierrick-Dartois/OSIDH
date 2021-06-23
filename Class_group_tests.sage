import sys
sys.path.append('Documents/Codes/OSIDH')
from Class_group import *

I=IdealClass(2,2,1)
J=IdealClass(1,0,1)
K=I*J
print(str(K))

L=IdealClass(-4)
print(str(L))
