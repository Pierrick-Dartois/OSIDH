import sys
sys.path.append('Documents/Codes/OSIDH')
from Class_group import *

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

## Ideal above a prime: dsic=-224=-2^2*56 (Disc_K=-56, conductor f=2)
# no 

# Test 7: 
