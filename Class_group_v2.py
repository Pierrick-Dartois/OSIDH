from sage.all import *

class IdealClassGroup:
	D_K_error="First argument D_K is not a fundamental discriminant"
	D_error="Argument is not a discriminant"

	def __init__(self,*args):
		# input of the form D_K (fundamental discriminant), f (conductor)
		if len(args)==2:
			D_K,f=args
			if D_K%4 not in [0,1]:
				raise ValueError(D_K_error)
			elif D_K%4==1:
				L=list(factor(D_K))
				for x in L:
					if x[1]>1:
						raise ValueError(D_K_error)
			else:
				m=D_K//4
				if m%4 not in [2,3]:
					raise ValueError(D_K_error)
				L=list(factor(m))
				for x in L:
					if x[1]>1:
						raise ValueError(D_K_error)
			self.D_K=D_K
			self.f=f
			self.D=D_K*f*f
			if D_K%4==0:
				K.<alpha>=NumberField(x**2-D_K//4)
			else:
				K.<alpha>=NumberField(x**2-x+(1-D_K)//4)
			self.K=K
		elif len(args)==1:
			[D]=args
			if D%4 not in [0,1]:
				raise ValueError(D_error)
			L=list(factor(D))
			D_K=1
			f=1
			for x in L:
				q,r=divmod(x[1],2)
				D_K*=x[0]**r
				f*=x[0]**q
			if D_K%4 not in [0,1]:
				D_K*=4
				f=f//2
			self.D_K=D_K
			self.f=f
			self.D=D
			if D_K%4==0:
				K.<alpha>=NumberField(x**2-D_K//4)
			else:
				K.<alpha>=NumberField(x**2-x+(1-D_K)//4)
			self.K=K

	def __str__(self):
		return "Ideal Class Group of discriminant {0}, fundamental discriminant {1} and conductor {2}".format(self.D,self.D_K,self.f)



class IdealClass:
	def __init__(self,*args):
		pass