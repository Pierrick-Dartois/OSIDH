from sage.quadratic_forms.binary_qf import BinaryQF
from sage.rings.all import Integers, ZZ
from sage.arith.all import gcd
from sage.functions.other import ceil
from sage.rings.finite_rings.integer_mod import Mod

class IdealClass:
	def __init__(self,*args):
		if len(args)==3:# Usual definition
			a,b,c=args
			d=b*b-4*a*c
			if a<0 or d>=0 or gcd([a,b,c])!=1:
				raise ValueError('This is not a positive definite primitive quadratic form')
			else:
				self.a=a
				self.b=b
				self.c=c
				self.disc=d
				self.reduce()
		elif len(args)==2:# Returns a prime ideal lying above a prime
			d,p=args# p shoud be a prime but it is not tested because it is costly
			if d>=0 or d%4 not in [0,1]:
				raise ValueError('The discriminant should be negative and congruent to 0 or 1 mod 4')
			elif Mod(d,p).is_square():
				R=Integers(4*p)
				b=ZZ(R(d).sqrt())
				self.a=p
				self.b=b
				self.c=(b*b-d)//(4*p)
				self.disc=d
				self.reduce()
			elif d%4==0:
				self.a=1
				self.b=0
				self.c=-d//4
				self.disc=d
			else:
				self.a=1
				self.b=1
				self.c=(1-d)//4
				self.disc=d
		elif len(args)==1:# Identity element
			[d]=args
			if d>=0:
				raise ValueError('The discriminant should be negative')
			if d%4==0:
				self.a=1
				self.b=0
				self.c=-d//4
				self.disc=d
			elif d%4==1:
				self.a=1
				self.b=1
				self.c=(1-d)//4
				self.disc=d
			else:
				raise ValueError('The discriminant should be congruent to 0 or 1 mod 4')

	def reduce(self):
		a,b,c=self.a,self.b,self.c
		while not -a<b<=a<=c or (b<0 and a==c):
			# Make sure a<=c
			if a>=c:
				a,c=c,a
				b=-b
			# Make sure abs(b)<= a
			if c>=a:
				if b>a:
					m=-ceil(b/(2*a))
					c=a*m*m+b*m+c
					b=b+2*a*m
				elif -b>a:
					m=ceil(-b/(2*a))
					c=a*m*m+b*m+c
					b=b+2*a*m
			# Make sure that b>=0 when |b|=a or a=c
			if (c>=a and b==-a) or (c==a and b<0):
				b=-b
		self.a,self.b,self.c=a,b,c

	def __str__(self):
		return "{0}*x^2 + {1}*x*y + {2}*y^2".format(self.a,self.b,self.c)

	def __eq__(self,q):
		return self.a==q.a and self.b==q.b and self.c==q.c

	def __ne__(self,q):
		return not self.__eq__(q)

	def __mul__(self,q):
		if self.disc!=q.disc:
			raise ValueError('The ideals are in distinct orders')
		else:
			a1=q.a
			b1=q.b
			c1=q.c
			if gcd([self.a,a1,(self.b+b1)//2])!=1:
				a1,c1=c1,a1
				b1=-b1
				if gcd([self.a,a1,(self.b+b1)//2])!=1:
					c1=a1+b1+c1
					b1=b1+2*a1
					a1,c1=c1,a1
					b1=-b1
			R=Integers(4*self.a*a1)
			L_B=R(self.disc).sqrt(all=True)
			L_B=[ZZ(x) for x in L_B]
			b_continue=True
			i=0
			while b_continue:
				if (L_B[i]-self.b)%(2*self.a)==0 and (L_B[i]-b1)%(2*a1)==0:
					B=L_B[i]
					b_continue=False
				else:
					i+=1
			return IdealClass(self.a*a1,B,(B*B-self.disc)//(4*self.a*a1))

	def inverse(self):
		return IdealClass(self.a,-self.b,self.c)

	def __pow__(self,k):
		l=max(-k,k)
		L=l.digits(2)
		n=len(L)
		q=IdealClass(self.disc)
		for i in range(n):
			if L[n-1-i]==0:
				q=q*q
			else:
				q=(q*q)*self
		if k>=0:
			return q
		else:
			return q.inverse()


	
