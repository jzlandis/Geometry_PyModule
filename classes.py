#!/usr/bin/python3

#defining classes that behave like vectors

from copy import deepcopy
import numpy as np
import math

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
			yield i
	def __getitem__(self,index):
		return self.E[index]

class _vector():
	def __init__(self,components,basis=None):
		'''components should be in system (defaults to None)'''
		if isinstance(components,self.__class__):
			if basis == None:
				self.comps = [float(x) for x in components._toGlobal()]
				self.basis = None
			else:
				self.comps = [float(x) for x in components]
				self.basis = deepcopy(basis)
		elif isinstance(components,list) or isinstance(components,tuple):
			self.comps = [float(x) for x in components]
			self.basis = deepcopy(basis)
		else:
			_x = components._toOtherType(self.__class__.__name__)
			self.comps = _x.comps
			self.basis = _x.basis
	def __iter__(self):
		for i in self.comps:
			yield i
	def __repr__(self):
		return "%s(%s)\n basis = %s" % (self.__class__.__name__,str(self.comps),str(self.basis))
	def __getitem__(self,index):
		return self.comps[index]

class sphericalVector(_vector):
	def _toOtherType(self,t):
		if t == "cartesianVector":
			r,theta,phi = self.comps
			x = r*math.cos(theta)*math.sin(phi)
			y = r*math.sin(theta)*math.sin(phi)
			z = r*math.cos(phi)
			return cartesianVector([x,y,z],basis=self.basis)
		else:
			print "Type changer didn't catch an option, was passed %s" % str(t)

class cartesianVector(_vector):
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
	def _toOtherType(self,t):
		if t == "sphericalVector":
			r = sum([x**2 for x in self.comps])**0.5
			theta = math.atan(self.comps[1]/self.comps[0])
			phi = math.acos(self.comps[2]/r)
			return sphericalVector([r,theta,phi],basis=self.basis)
		else:
			print "Type changer didn't catch an option, was passed %s" % str(t)
	def _toGlobal(self):
		return self._changeBasis(GLOBAL_CARTESIAN_BASIS)
	def magnitude(self):
		return sum([x**2 for x in self.comps])**0.5
	def unitize(self):
		return cartesianVector([x / self.magnitude() for x in self],basis=self.basis)
	def dot(self,other):
		#return sum([x*y for x,y in zip(self._toGlobal(),cartesianVector(other)._toGlobal())])
		return sum([x*y for x,y in zip(self,other._changeBasis(self.basis))])
	def cross(self,other):
		otherChangedBasis = other._changeBasis(self.basis)
		return cartesianVector([self[i]*otherChangedBasis[j]-self[j]*otherChangedBasis[i] for i,j in ((1,2),(2,0),(0,1))],basis=self.basis)
	def angle(self,other):
		return math.acos(self.dot(other._changeBasis(self.basis)) / (self.magnitude()*other.magnitude()))
	def area(self,other):
		return 0.5*self.cross(other).magnitude()
	def volume(self,other1,other2):
		return abs(other2.dot(self.cross(other1)))/6.
	def normal(self,other):
		return cartesianVector(cartesianBasis((self,other))[2])._changeBasis(self.basis)
	def __add__(self,other):
		return cartesianVector([x+y for x,y in zip(self,other._changeBasis(self.basis))],basis=self.basis)
	def __sub__(self,other):
		return cartesianVector([x-y for x,y in zip(self,other._changeBasis(self.basis))],basis=self.basis)

GLOBAL_CARTESIAN_BASIS = cartesianBasis(((1.,0.,0.),(0.,1.,0.)))

class cartesianScalarField():
	def __init__(self,funcs,basis=None,origin=(0.,0.,0.)):
		self.funcs = funcs
	def __call__(self,position):
		return sum([x(position[i]) for i,x in enumerate(self.funcs)])


yAxisBasis = cartesianBasis([[0.,1.,0.],[0.,0.,1.]])
yAxisVector = cartesianVector([1,0,0],basis=yAxisBasis)
xAxisBasis = cartesianBasis([[1,0,0],[0,1,0]])
xAxisVector = cartesianVector([1,0,0],basis=xAxisBasis)
zAxisBasis = cartesianBasis([[0.,0.,1.],[1.,0.,0.]])
zAxisVector = cartesianVector([1,0,0],basis=zAxisBasis)

def funca(x):
	return x**0.2
def funcb(x):
	return x/4.
def funcc(x):
	return x+2

print("Should be a global z:")
print(xAxisVector.cross(yAxisVector))
print("Should be zero:")
print(xAxisVector.dot(yAxisVector))
print("Should be 90 degrees (1.57 rad):")
print(xAxisVector.angle(yAxisVector))
print("Should be 0.5 (triangular area between yAxisVector and yAxisVector)")
print(xAxisVector.area(yAxisVector))
print("Should be 1/6 (0.1667) (tetrahedron volume between x,y,zAxisVectors)")
print(xAxisVector.volume(yAxisVector,zAxisVector))
print("Should be a global z ")
print(yAxisVector.normal(xAxisVector)._toGlobal())

cSFex = cartesianScalarField((funca,funcb,funcc))
print(cSFex((1,1,1)))

print yAxisVector
print sphericalVector(yAxisVector)
print cartesianVector(sphericalVector(yAxisVector))
