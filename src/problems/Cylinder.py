__author__ = "Kristian Valen-Sendstad <kvs@simula.no>"
__date__ = "2009-10-01"
__copyright__ = "Copyright (C) 2009-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg, 2010.
# Modified by Erich L Foster, 2013

from problembase import *
from numpy import array

# Constants related to the geometry
bmarg = 1.e-3 + DOLFIN_EPS
xmin = 0.0
xmax = 2.2
ymin = 0.0
ymax = 0.41
xcenter = 0.2
ycenter = 0.2
radius = 0.05

class InitialConditions(Expression):
    def eval(self,values,x):
        values[0] = 0.
        values[1] = 0.
        values[2] = 0.

    def value_shape(self):
      return (3,)

# Inflow boundary
class InflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] < xmin + bmarg

# No-slip boundary
class NoslipBoundary(SubDomain):
    def inside(self, x, on_boundary):
        dx = x[0] - xcenter
        dy = x[1] - ycenter
        r = sqrt(dx*dx + dy*dy)
        return on_boundary and \
               (x[1] < ymin + bmarg or x[1] > ymax - bmarg or \
                r < radius + bmarg)

# Outflow boundary
class OutflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] > xmax - bmarg

# Problem definition
class Problem(ProblemBase):

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Load mesh
        Nx = options["Nx"]
        Ny = options["Ny"]

        rect = Rectangle(xmin, ymin, xmax, ymax)
        circ = Circle(xcenter, ycenter, radius)
        domain = rect - circ
        self.mesh = Mesh(domain, Nx)

        self.t0 = 0.
        self.T = options['T']
        self.k = options['dt']

    def initial_conditions(self, W):
        w0 = InitialConditions()
        w0 = project(w0,W)

        return w0

    def boundary_conditions(self, W, t):
        # Create inflow boundary condition
        g0 = Expression(('4*Um*(x[1]*(ymax-x[1]))/(ymax*ymax)*4*t*t/(4*t*t+1)', '0.0'), Um=1.5, ymax=ymax, t=t)

        bc0 = DirichletBC(W.sub(0), g0, InflowBoundary())

        # Create no-slip boundary condition
        bc1 = DirichletBC(W.sub(0), Constant((0,0)), NoslipBoundary())

        # Create outflow boundary condition for pressure
        bc2 = DirichletBC(W.sub(1), Constant(0), OutflowBoundary())

        # Collect boundary conditions
        bcs = [bc0, bc1, bc2]

        return bcs

    def F1(self, t):
        #forcing function for the momentum equation
        return Constant((0,0))

    def F2(self, t):
        #mass source for the continuity equation
        return Constant(0)

    def __str__(self):
        return "Cylinder"
