__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"
#
#   adapted from channel.py in nsbench originally developed 
#   by Kent-Andre Mardal <kent-and@simula.no>
#

'''
This is an extremely boring problem with no forcing and on a square.
This is basically a blank problem that we can adapt with optional inputs.
'''

from problembase import *
from numpy import array

class InitialConditions(Expression):
    def __init__(self,p1,p2):
        self.p1 = p1
        self.p2 = p2
        self.Ap = 1.
        self.An = 1.
        self.s = 1E-8

    def eval(self,values,x):
        xp = self.p1[0]; yp = self.p1[1]; 
        xn = self.p2[0]; yn=self.p2[1] 
        Ap = self.Ap; An = self.An
        s = self.s

        values[0] = 0.
        values[1] = 0.
        values[2] = Ap*exp(-((x[0]-xp)**2. + (x[1]-yp)**2.)/(2.*s*s)) \
            + An*exp(-((x[0] - xn)**2. + (x[1] - yn)**2.)/(2.*s*s))

    def value_shape(self):
      return (3,)

# No-slip boundary
class NoslipBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary 

# Problem definition
class Problem(ProblemBase):
#   2D channel flow.

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Create mesh
        Nx = options["Nx"]
        Ny = options["Ny"]
        x0 = 0
        x1 = 1
        y0 = 0
        y1 = 1
        self.p1 = Point(x0 + 0.25*(x1-x0),y0 + 0.5*(y1-y0))
        self.p2 = Point(x1 - 0.25*(x1-x0),y0 + 0.5*(y1-y0))
        self.mesh = RectangleMesh(x0,y0,x1,y1,Nx, Ny)

        self.t0 = 0.
        self.T = options['T']
        self.k = options['dt']

    def initial_conditions(self, W):
        w0 = InitialConditions(self.p1, self.p2)
        w0 = project(w0,W)

        return w0

    def boundary_conditions(self, W, t):
        # Create no-slip boundary condition for velocity
        bcs = DirichletBC(W.sub(0), Constant((0.0, 0.0)), NoslipBoundary())

        return bcs

    def F1(self, t):
        #forcing function for the momentum equation
        return Expression(self.options['F1'],t=t)

    def F2(self, t):
        #mass source for the continuity equation
        #our point source is defined in the initial condition since I need the
        #functional space to apply it
        return Expression('Ap*exp(-(pow(x[0]-xp,2) + pow(x[1]-yp, 2))/(2*s*s)) + An*exp(-(pow(x[0] - xn,2) + pow(x[1] - yn, 2))/(2*s*s))', xp=self.p1[0], yp=self.p1[1], xn=self.p2[0], yn=self.p2[1], Ap=1.0, An=-1.0, s=1E-8, t=t)

    def __str__(self):
        return 'Point'
