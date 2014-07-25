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

bmarg = 1.e-3 + DOLFIN_EPS
xmin = -1.
xmax = 1.
ymin = -1.
ymax = 1
xcenter = 0.
ycenter = 0.
radius = 0.3

A = 10.
kappa = 100;
rho = 1500.; c = 1480

class InitialConditions(Expression):
    def eval(self,values,x):
        values[0] = 0.

class OuterBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (near(x[0],xmin) or near(x[0],xmax) \
                or near(x[1],ymin) or near(x[1],ymax))

class InnerBoundary(SubDomain):
    def inside(self, x, on_boundary):
        dx = x[0] - xcenter
        dy = x[1] - ycenter
        r = sqrt(dx*dx + dy*dy)
        return on_boundary and r < radius + bmarg


# Problem definition
class Problem(ProblemBase):
#   2D channel flow.

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Create mesh
        self.Nx = options["Nx"]

        rect = Rectangle(xmin, ymin, xmax, ymax)
        circ = Circle(xcenter, ycenter, radius)
        domain = rect - circ
        self.mesh = Mesh(domain, self.Nx)

        self.t0 = 0.
        self.T = options['T']
        self.k = options['dt']

        self.kappa = kappa
        self.rho = rho
        self.c = c

    def initial_conditions(self, W):
        w0 = InitialConditions()
        w0 = project(w0,W)

        return w0

    def boundary_conditions(self, W, t):
        # Create no-slip boundary condition for velocity
        bc0 = DirichletBC(W, Constant(0.0), OuterBoundary())
        bc1 = DirichletBC(W, A, InnerBoundary())

        return [bc0, bc1]

    def F1(self, t):
        #forcing function for the momentum equation
        return Expression(self.options['F1'],t=t)

    def F2(self, t):
        #mass source for the continuity equation
        return Expression(self.options['F2'],t=t)

    def __str__(self):
        return 'Plate'
