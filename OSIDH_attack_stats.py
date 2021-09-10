import sys
sys.path.append("Documents/Codes/OSIDH")
from sage.all import *
import csv
from OSIDH_protocol import *
from OSIDH_attack import *

# Parameters
n=28
t=10
l=2
r=3
d_K=-4
d=d_K*l**(2*n)
N_sample=10

osidh=OSIDH(n,t,l,r,d_K)
pub_chain=Chain(osidh)
L_keys_A,L_keys_B,L_protocol,L_step1a,L_step1b,L_step2,L_step3,L_attack=OSIDH_attack_stats(osidh,pub_chain,N_sample)

with open("Documents/Codes/OSIDH/Data_files/Attack_stats.csv","w",encoding="utf-8") as f:
	writer=csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
	writer.writerow(["Key Alice"]+[""]*(t-1)+["Key Bob"]+[""]*(t-1)+["Time protocol","Time step 1 - Alice","Time step 1 - Bob","Time step 2","Time step 3","Total attack time"])
	for i in range(N_sample):
		writer.writerow(L_keys_A[i]+L_keys_B[i]+[L_protocol[i],L_step1a[i],L_step1b[i],L_step2[i],L_step3[i],L_attack[i]])
