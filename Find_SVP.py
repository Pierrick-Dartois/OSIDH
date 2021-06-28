import sys
sys.path.append('Documents/Codes/OSIDH')
from sage.all import *
from fpylll import *

with open("Documents/Codes/OSIDH/Relations_lattice_basis.txt","r",encoding="utf-8") as f:
	L=f.readlines()
	i0=1
	while L[i0]!="Lattice basis\n":
		i0+=1
	i0+=1
	t=int(L[i0][0:-1])
	i0+=1
	M=[]
	for i in range(t):
		x=L[i+i0][:-1]
		row=x.split(",")
		M.append([int(y) for y in row])

M=IntegerMatrix.from_matrix(M).transpose()
M=LLL.reduction(M)
u=SVP.shortest_vector(M)
print(u)
#M_bkz = BKZ.reduction(M.transpose(), BKZ.Param(3))
#L_max=[]
#for i in range(74):
	#Li=list(M_bkz[i])
	#L_max.append(max(max(Li),-min(Li)))
#print(min(L_max))
