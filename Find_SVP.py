import sys
sys.path.append('Documents/Codes/OSIDH')
from sage.all import *
from fpylll import *
from time import time

with open("Documents/Codes/OSIDH/Data_files/Relations_lattice_basis.txt","r",encoding="utf-8") as f:
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

b_LLL=False
k_BKZ=4

if b_LLL:
	t1=time()
	M_red=LLL.reduction(M)
	t2=time()
	print("Time LLL: {0} s".format(t2-t1))
else:
	t1=time()
	M_red = BKZ.reduction(M, BKZ.Param(k_BKZ))
	t2=time()
	print("Time BKZ: {0} s".format(t2-t1))

# Displaying the shortest vector
L_max=[]
for i in range(t):
	L_max.append(max(max(M_red[i]),-min(M_red[i])))
N_inf=min(L_max)
i0=L_max.index(N_inf)
u=M_red[i0]
print("Infinity norm of the shortest vector: {0}".format(N_inf))
print("Shortest vector:")
print(u)

# Saving the LLL/BKZ matrix (outputted lattice vectors are rows, unlike entry lattice vectors of M given in columns) 
with open("Documents/Codes/OSIDH/Data_files/Lattice_BKZ_4.txt","w",encoding="utf-8") as f:
	f.write("Infinity norm of the shortest vector:\n")
	f.write(str(N_inf)+"\n")
	f.write("Index of the shortest vector:\n")
	f.write(str(i0)+"\n")
	f.write("Shortest vector:\n")
	f.write(str(M_red[i0])[1:-1]+"\n")
	f.write("Reduced matrix:\n")
	f.write(str(t)+"\n")
	for i in range(t):
		f.write(str(M_red[i])[1:-1]+"\n")
