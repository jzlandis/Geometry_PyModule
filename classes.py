#!/usr/bin/python3

#defining classes that behave like vectors

from copy import deepcopy
import numpy as np

def dot(v1,v2):
	return sum([x*y for x,y in zip(v1,v2)])

class cartesianBasis():
	def __init__(self,basisVectors):
		#self.E = basisUnitVectors
		if len(basisVectors) == 2:
			temp_e1,temp_e2 = [cartesianVector(y).unitize() for y in basisVectors]
			temp_e3 = temp_e1.cross(temp_e2)
			self.E = [list(temp_e1),list(temp_e2),list(temp_e3)]
			del temp_e1,temp_e2,temp_e3
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
			q = np.array([[cartesianVector(newBasis[i]).dot(cartesianVector(refBasis[j])) for i in (0,1,2)] for j in (0,1,2)])
			a = np.array([[x,] for x in list(self)])
			return cartesianVector([sum(x) for x in (a*q).T],basis=newBasis)
	def _toGlobal(self):
		return self._changeBasis(GLOBAL_CARTESIAN_BASIS)
	def magnitude(self):
		return sum([x**2 for x in self.comps])**0.5
	def unitize(self):
		return cartesianVector([x / self.magnitude() for x in self],basis=self.basis)
	def __repr__(self):
		return "cartesianVector(%s)\n basis = %s" % (str(self.comps),str(self.basis))
	def __iter__(self):
		for i in self.comps:
			yield i
	def dot(self,other):
		#return sum([x*y for x,y in zip(self._toGlobal(),cartesianVector(other)._toGlobal())])
		return sum([x*y for x,y in zip(self,other._changeBasis(self.basis))])
	def cross(self,other):
		otherChangedBasis = other._changeBasis(self.basis)
		return cartesianVector([self[i]*otherChangedBasis[j]-self[j]*otherChangedBasis[i] for i,j in ((1,2),(2,0),(0,1))],basis=self.basis)
	def __add__(self,other):
		return cartesianVector([x+y for x,y in zip(self,other._changeBasis(self.basis))],basis=self.basis)
	def __sub__(self,other):
		return cartesianVector([x-y for x,y in zip(self,other._changeBasis(self.basis))],basis=self.basis)
	def __getitem__(self,index):
		return self.comps[index]
GLOBAL_CARTESIAN_BASIS = cartesianBasis(((1.,0.,0.),(0.,1.,0.)))

yAxisBasis = cartesianBasis([[0.,1.,0.],[0.,0.,1.]])
yAxisVector = cartesianVector([1,0,0],basis=yAxisBasis)
xAxisBasis = cartesianBasis([[1,0,0],[0,1,0]])
xAxisVector = cartesianVector([1,0,0],basis=xAxisBasis)

print "Should be a global z:"
print xAxisVector.cross(yAxisVector)
print "Should be zero:"
print xAxisVector.dot(yAxisVector)
