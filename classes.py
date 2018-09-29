#!/usr/bin/python3

#defining classes that behave like vectors

from copy import deepcopy
import numpy as np

def dot(v1,v2):
	return sum([x*y for x,y in zip(v1,v2)])

class cartesianBasis():
	def __init__(self,basisUnitVectors):
		#self.E = basisUnitVectors
		if len(basisUnitVectors) == 3:
			self.E = [list(cartesianVector(y).unitize()) for y in basisUnitVectors]
		#check E for linear independence
	def __repr__(self):
		return "cartesianBasis(%s)" % ", ".join([str("e%d = " % (i+1)) + str(self.E[i]) for i in (0,1,2)])
	def __eq__(self,other):
		if isinstance(other,cartesianBasis):
			return all([all([x==y for x,y in zip(s,t)]) for s,t in zip(other.E,self.E)])
		else:
			return False
	def __iter__(self):
		for i in self.E:
			return i
	def __getitem__(self,index):
		return self.E[index]

def cartBasisFromPoints(points):
	'''direction from first to second point defines e1 in resulting basis, third point defines a point on the e1-e2 plane'''
	p1,p2,p3 = [[float(x) for x in y] for y in points]
	e1 = [y-x for x,y in zip(p1,p2)]

class cartesianVector():
	def __init__(self,components,basis=None):
		self.comps = [float(x) for x in components]
		self.basis = deepcopy(basis)
	def _changeBasis(self,newBasis):
		if newBasis == self.basis:
			return self
		else:
			if self.basis == None:
				refBasis = GLOBAL_CARTESIAN_BASIS
			else:
				refBasis = self.basis
			q = np.array([[cartesianVector(newBasis[j]).dot(cartesianVector(refBasis[i])) for i in (0,1,2)] for j in (0,1,2)])
			a = np.array([[x,] for x in list(self)])
			return cartesianVector([sum(x) for x in a*q],basis=newBasis)
	def _toGlobal(self):
		return self._changeBasis(GLOBAL_CARTESIAN_BASIS)
	def magnitude(self):
		return sum([x**2 for x in self.comps])**0.5
	def unitize(self):
		return cartesianVector([x / self.magnitude() for x in self])
	def __repr__(self):
		return "cartesianVector(%s)\n basis = %s" % (str(self.comps),str(self.basis))
	def __iter__(self):
		for i in self.comps:
			yield i
	def dot(self,other):
		#return sum([x*y for x,y in zip(self._toGlobal(),cartesianVector(other)._toGlobal())])
		return sum([x*y for x,y in zip(self,other._changeBasis(self.basis))])
	def __add__(self,other):
		return cartesianVector([x+y for x,y in zip(self,other._changeBasis(self.basis))],basis=self.basis)
	def __sub__(self,other):
		return cartesianVector([x-y for x,y in zip(self,other._changeBasis(self.basis))],basis=self.basis)
GLOBAL_CARTESIAN_BASIS = cartesianBasis(((1,0,0),(0,1,0),(0,0,1)))
b = cartesianVector([1,2,3])
b1 = cartesianVector([1,2,3])
b2 = cartesianVector([2,3,4])
c = cartesianBasis(((0.5,3,2),(4,5,3),(1,2,3)))
c1 = cartesianVector([1,2,3],basis=c)
ls = dir()
print b._toGlobal()
print b._changeBasis(c)
