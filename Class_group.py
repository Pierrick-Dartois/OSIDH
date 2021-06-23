from sage.quadratic_forms.binary_qf import BinaryQF
from sage.rings.all import Integers, ZZ

class IdealClass:
	def __init__(self,*args):
		if len(args)==3:
			a,b,c=args
			self.a=a
			self.b=b
			self.c=c
			self.disc=b*b-4*a*c
		elif len(args)==1:
			[d]=args
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

	def __str__(self):
		return "{0}*x^2 + {1}*x*y + {2}*y^2".format(self.a,self.b,self.c)

	def __mul__(self,q):
		if self.disc!=q.disc:
			raise ValueError('The ideals are in distinct orders')
		else:
			R=Integers(4*self.a*q.a)
			L_B=R(self.disc).sqrt(all=True)
			L_B=[ZZ(x) for x in L_B]
			b_continue=True
			i=0
			while b_continue:
				if (L_B[i]-self.b)%(2*self.a)==0 and (L_B[i]-q.b)%(2*q.a)==0:
					B=L_B[i]
					b_continue=False
				else:
					i+=1
			return IdealClass(self.a*q.a,B,(B*B-self.disc)//(4*self.a*q.a))


	
