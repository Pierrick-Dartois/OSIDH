import sys
from Class_group import *
from sage.all import *

### Testing the constructor
print("Testing the constructor")

## Definition with 3 arguments

# Test 1: two forms of disc -4
try:
	I=IdealClass(2,2,1)
	J=IdealClass(1,0,1)
	print("Test 1: success")
except:
	print("Test 1:failure")

# Test 2: a form of discriminant 12
try:
	I=IdealClass(1,4,1) 
	print("Test 2: failure")
except ValueError:
	print("Test 2: success")

## Identity element

# Test 3: identity element with disc congruent to 0 mod 4
try:
	E=IdealClass(-24)
	print("Test 3: success")
except:
	print("Test 3: failure")

# Test 4: identity element with disc congruent to 1 mod 4
try:
	E=IdealClass(-3)
	print("Test 4: success")
except:
	print("Test 4: failure")

# Test 5: identity element with positive discriminant
try:
	E=IdealClass(4)
	print("Test 5: failure")
except ValueError:
	print("Test 5: success")

# Test 6: identity element with discriminant not congruent to 0 or 1 mod 4
try:
	E=IdealClass(-25)
	print("Test 6: failure")
except ValueError:
	print("Test 6: success")

## Ideal above a prime: dsic=-224=-2**2*56 (Disc_K=-56, conductor f=2)

# Test 7: splitting prime 5 ((Delta_K/5)=+1)
try:
	P=IdealClass(-224,5)
	print("Test 7: success. Creation of: P="+str(P)+" lying above 5 (split for D=-224) complete.")
except:
	pint("Test 7: failure")

# Test 8: ramifying prime 7 ((Delta_K/7)=0)
try:
	Q=IdealClass(-224,7)
	print("Test 8: success. Creation of: Q="+str(Q)+" lying above 7 (ramifying for D=-224) complete.")
except:
	pint("Test 8: failure")

# Test 9: inert prime 11 ((Delta_K/11)=-1)
try:
	R=IdealClass(-224,11)
	print("Test 9: success. Creation of: R="+str(R)+" lying above 11 (inert for D=-224) complete.")
except:
	pint("Test 9: failure")

# Test 10: inert prime 5 for Delta_K= 1 mod 4 ((Delta_K/5)=-1)
try:
	R=IdealClass(-48,5)
	print("Test 10: success. Creation of: R="+str(R)+" lying above 5 (inert for D=-48) complete.")
except:
	pint("Test 10: failure")

### Testing the accuracy of the reduction method (applied systematically when an object is created)
print("\nTesting reduction method")

# Test 1: usual definition
b_test=True
for i in range(1000):
	a=randint(1,1000)
	b=randint(-1000,1000)
	c=randint(1,1000)
	if b*b-4*a*c<0 and gcd([a,b,c])==1:
		I=IdealClass(a,b,c)
		a1,b1,c1=I.a,I.b,I.c
		q=BinaryQF(a1,b1,c1)
		if not q.is_reduced():
			b_test=False

if b_test:
	print("Test 1: success")
else:
	print("Test 1: failure")

# Test 2: identity
b_test=True
for i in range(1000):
	d=randint(-1000,-3)
	if d%4 in [0,1]:
		I=IdealClass(d)
		a1,b1,c1=I.a,I.b,I.c
		q=BinaryQF(a1,b1,c1)
		if not q.is_reduced():
			b_test=False

if b_test:
	print("Test 2: success")
else:
	print("Test 2: failure")

# Test 3: prime ideals
b_test=True
for i in range(1000):
	d=randint(-1000,-3)
	p=random_prime(1000)
	if d%4 in [0,1] and d%p!=0:
		I=IdealClass(d,p)
		a1,b1,c1=I.a,I.b,I.c
		q=BinaryQF(a1,b1,c1)
		if not q.is_reduced():
			b_test=False

if b_test:
	print("Test 3: success")
else:
	print("Test 3: failure")

### Testing equality
print("\nTesting == and != operators")

# Test 1: R of test 10 is identity
if R==IdealClass(-48):
	print("Test 1: success")
else:
	print("Test 1: failure")

# Test 2: Q of test 8 is not identity
if Q!=IdealClass(-224):
	print("Test 2: success")
else:
	print("Test 2: failure")

### Testing multiplication
print("\nTesting multiplication")

# Test 1: it should work (same discriminants)
d=-4*3**10
L_forms=BinaryQF_reduced_representatives(d,primitive_only=True)
n=len(L_forms)
b_test=True
for i in range(100):
	j=randint(0,n-1)
	k=randint(0,n-1)
	f=L_forms[j]
	g=L_forms[k]
	I=IdealClass(f[0],f[1],f[2])
	J=IdealClass(g[0],g[1],g[2])
	try:
		K=I*J
	except:
		b_test=False

if b_test:
	print("Test 1: success")
else:
	print("Test 1: failure")

# Test 2: it should not work (distinct discriminants)
I=IdealClass(-56,5)
j=randint(0,n-1)
f=L_forms[j]
J=IdealClass(f[0],f[1],f[2])

try:
	K=I*J
	print("Test 2: failure")
except:
	print("Test 2: success")

## Testing exponentiation
print("\nTesting exponentiation")

# Test 1: recovering neutral element
E=IdealClass(d)
if J**n==E:
	print("Test 1: success")
else:
	print("Test 1: failure")

# Test 2: huge exponentiation to check special cases and timing
d=-3*2**100
I=IdealClass(d,7)
k=randint(1,2**100)
l=-randint(1,2**100)
try:
	J=I**k
	K=I**l
	print("Test 2: success")
except:
	print("Test 2: failure")
