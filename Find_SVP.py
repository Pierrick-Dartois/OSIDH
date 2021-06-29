import sys
sys.path.append('Documents/Codes/OSIDH')
from sage.all import *
from fpylll import *
from time import time

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

M=IntegerMatrix.from_matrix(M)

t1=time()
M_bkz = BKZ.reduction(M.transpose(), BKZ.Param(42))
t2=time()
print("Time BKZ: {0} s".format(t2-t1))

# Saving the BKZ matrix (outputted lattice vectors are rows, unlike entry lattice vectors of M given in columns) 
with open("Documents/Codes/OSIDH/Lattice_BKZ_42.txt","w",encoding="utf-8") as f:
	f.write(str(t)+"\n")
	for i in range(t):
		f.write(str(M_bkz[i])[1:-1]+"\n")