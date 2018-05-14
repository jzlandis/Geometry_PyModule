#!/usr/bin/python
import sys
from itertools import combinations

MODULE_NAME = "landis_geometry"
MODULE_ERROR = "Error in %s module:" % MODULE_NAME

class vector:
    '''landis_geometry vector class'''
    def __init__(self,ins,orig = None,typ = 0,syst = 0):
        if not type(ins) == list or sum([1 for x in ins if not type(x) in [float,int]]) > 0:
            sys.stderr.write("%s vector class was initialized with the following input where a list of floats or ints should be: %s \"%s\"\n" % (MODULE_ERROR,str(type(ins)),ins))
            exit()
        self.coords = ins
        if typ == 0:
            self.kind = "cartesian"
        else:
            self.kind = typ
        if not orig == None:
            if not type(orig) == list or sum([1 for x in orig if not type(x) in [float,int]]) > 0:
                sys.stderr.write("%s vector class was initialized with the following input where a list of floats or ints should be: %s \"%s\"\n" % (MODULE_ERROR,str(type(orig)),orig))
                exit()
            else:
                self.origin = orig
        else:
            self.origin = None
        if syst == 0:
            self.system = "global"
        else:
            self.system = syst
    def __repr__(self):
        return "vector(%s, origin=%s, typ=\"%s\", syst=\"%s\")" % (self.coords,self.origin,self.kind,self.system)
    def __iter__(self):
        for j in self.coords:
            yield j
    def length(self):
        '''returns length of the vector'''
        return sum([x**2 for x in self.coords]) ** 0.5
    def __len__(self):
        '''returns dimension of the vector'''
        return len(self.coords)
    def __getitem__(self,i):
        '''returns indexed coord'''
        return self.coords[i]
    def unit(self):
        '''returns unit vector version of the vector'''
        return vector([x/self.length() for x in self.coords],orig=self.origin,syst=self.system)
    def terminus(self):
        '''returns the point that the vector points to if the vector has an origin'''
        if not self.origin == None:
            return [x+y for x,y in zip(self.origin,self.coords)]
        else:
            return None
    def to_system(self,system,orig_move=True):
        if self.system == system:
            return self
        else:
            if system.dimension == 0:
                if orig_move:
                    return vector(self.coords,system.origin,syst=system)
                else:
                    return vector(self.coords,self.origin,syst=system)
            elif system.dimension == 1:
                if orig_move:
                    return vector([dot(self.coords,system.x_hat)],system.origin,syst=system)
                else:
                    return vector([dot(self.coords,system.x_hat)],self.origin,syst=system)
            elif system.dimension == 3:
                if orig_move:
                    return vector([dot(self.coords,i) for i in [system.x_hat,system.y_hat,system.z_hat]],system.origin,syst=system)
                else:
                    return vector([dot(self.coords,i) for i in [system.x_hat,system.y_hat,system.z_hat]],self.origin,syst=system)
    def move_to_system(self,system,orig_move=True):
        self.system = system
        if system.dimension >= 0:
            if orig_move:
                self.origin = system.origin
            else:
                pass
        if system.dimension == 1:
            self.coords = [dot(self.coords,system.x_hat)]
        elif system.dimension == 3:
            self.coords = [dot(self.coords,i) for i in [system.x_hat,system.y_hat,system.z_hat]]
    def to_global(self):
        if self.system == "global":
            return vector(self.coords)
        else:
            if self.origin == None:
                _o = self.system.origin
            else:
                _o = [a+b for a,b in zip(self.system.origin,p2cs(self.origin,self.system.vecs))]
            if self.system.dimension == 0:
                return vector(self.coords,orig=self.origin)
            elif self.system.dimension == 1:
                return vector([self.coords[0]*self.system.x_hat[x] for x in range(len(self.system.x_hat))],orig=self.system.origin)
            elif self.system.dimension == 3:
                #tempcoords = [0.0,0.0,0.0]
                #for i in range(len(self.coords)):
                    #if i == 0:
                #    tempcoords = [a+b for a,b in zip(tempcoords,[self.coords[i]*[self.system.x_hat[x],self.system.y_hat[x],self.system.z_hat[x]][i] for x in range(len(self.system.x_hat))])]
                #    #tempcoords.append([self.coords[i]*[self.system.x_hat[x],self.system.y_hat[x],self.system.z_hat[x]][i] for x in range(len(self.system.x_hat))])
                #tempcoords += [0.0,]*(3-len(tempcoords))
                #return vector(tempcoords,orig=self.system.origin)
                return vector(p2cs(self.coords,self.system.vecs),orig=self.system.origin)

def p2cs(p,cs):
    tempcoords = []
    for i in range(len(cs)):
        try:
            tempcoords.append([p[i]*cs[i][x] for x in range(len(cs[i]))])
        except IndexError:
            pass
        #tempcoords = [a+b for a,b in zip(tempcoords,[p[i]*cs[i][x] for x in range(len(cs[0]))])]
    return [sum([tempcoords[k][l] for k in range(len(tempcoords))]) for l in range(len(tempcoords[0]))]

def cross(v1,v2):
   i1,j1,k1 = v1
   i2,j2,k2 = v2
   i = j1*k2 - j2*k1
   j = -1 * ( i1*k2 - k1*i2 )
   k = i1*j2 - j1*i2
   del i1,j1,k1,i2,j2,k2
   return [i,j,k]
def dot(v1,v2):
   return sum([x*y for x,y in zip(v1,v2)])
def cart_colinear(*points):
    '''Determines whether all cartesian points input are colinear'''
    colinear = True
    if len(points) > 2:
        x,y,z = 0,1,2
        for point_combo in list(combinations(points,3)):
            x1,y1,z1 = point_combo[0]
            x2,y2,z2 = point_combo[1]
            x3,y3,z3 = point_combo[2]
            #print x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
            if x1*(y2*z3-z2*y3) - y1*(x2*z3 - z2*x3) + z1*(x2*y3-y2*x3) != 0:
            #if x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2) != 0:
                colinear = False
                break
    return colinear

class cart_system:
    def __init__(self,*sys_ins):
        for j in sys_ins:
            if not type(j) == list or sum([1 for x in j if not type(x) in [float,int]]) > 0:
                sys.stderr.write("%s system class was initialized with the following input where a list of floats or ints should be: %s \"%s\"\n" % (MODULE_ERROR,str(type(j)),j))
                exit()
        if len(sys_ins) >= 1:
            self.p1_in = sys_ins[0]
            self.origin = sys_ins[0]
            self.dimension = 0
        if len(sys_ins) >= 2:
            self.p2_in = sys_ins[1]
            if len(sys_ins[0]) == len(sys_ins[1]):
                self.x = vector([b-a for b,a in zip(sys_ins[1],self.origin)])
                self.x_hat = self.x.unit()
                self.dimension = 1
            else:
                sys.stderr.write("%s system class was initialized with the following points that are not of the same length: %s %s\n" % (MODULE_ERROR,sys_ins[0],sys_ins[1]))
                exit()
        if len(sys_ins) == 3:
            self.p3_in = sys_ins[2]
            if len(sys_ins[0]) == len(sys_ins[2]):
                if not cart_colinear(sys_ins[0],sys_ins[1],sys_ins[2]):
                    self.temp_y = vector([c-a for c,a in zip(sys_ins[2],sys_ins[0])])
                    self.z = vector(cross(self.x,self.temp_y))
                    self.z_hat = self.z.unit()
                    self.y_hat = vector(cross(self.z_hat,self.x_hat))
                    self.y = vector([self.temp_y.length()*x for x in self.y_hat])
                    del self.temp_y
                    self.dimension = 3
                else:
                    sys.stderr.write("%s system class was initialized with the following points that are colinear and thus cannot define a system: %s %s %s\n" % (MODULE_ERROR,sys_ins[0],sys_ins[1],sys_ins[2]))
                    #exit()
            else:
                sys.stderr.write("%s system class was initialized with the following points that are not of the same length: %s %s %s\n" % (MODULE_ERROR,sys_ins[0],sys_ins[1],sys_ins[2]))
                #exit()
        self.vecs = []
        if self.dimension >= 1:
            self.vecs.append(self.x_hat)
        if self.dimension >= 2:
            self.vecs.append(self.y_hat)
        if self.dimension == 3:
            self.vecs.append(self.z_hat)
    def __repr__(self):
        if self.dimension == 0:
            return "cart_system(%s) dimension = %s, origin at %s" % (self.origin,self.dimension,self.origin)
        elif self.dimension == 1:
            return "cart_system(%s,%s) dimension = %s, origin at %s" % (self.p1_in,self.p2_in,self.dimension,self.origin)
        elif self.dimension == 3:
            return "cart_system(%s,%s,%s) dimension = %s, origin at %s" % (self.p1_in,self.p2_in,self.p3_in,self.dimension,self.origin)
        else:
            return "ERROR"
    def __iter__(self):
        for j in self.vecs:
            yield j
    def __getitem__(self):
        return self.vecs[i]
    def __eq__(self,other):
        if self.__class__.__name__ == other.__class__.__name__:
            return [list(x) for x in self.vecs] + self.origin == [list(x) for x in other.vecs] + other.origin
        else:
            return False
a = [1,2,3]
b = [4,5,6]
c = [7,8,9]
c2 = [7,8,9.1]
v = vector([1,2,3])
s = cart_system(a,b,c2)
x = [1,0,0]
y = [0,1,0]
z = [0,0,1]
s2 = cart_system(x,y,z)
v2 = vector([1],syst=s2)
v2.to_global()
