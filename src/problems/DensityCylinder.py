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
xmin = 0.0; xmax = 2.2
ymin = 0.0; ymax = 0.41
xcenter = 0.2; ycenter = 0.2
radius = 0.05; Diameter = 2.*radius
Um = 1.5 #max velocity

At = 0.5 #Atwood number
rhoMin = 1.
rhoMax = rhoMin*(1. + At)/(1. - At)
S = 1./16.

class InitialConditions(Expression):
    def __init__(self):
        self.A1 = rhoMin
        self.A2 = rhoMax
        self.S = S #variance

    def eval(self,values,x):
        A1 = self.A1
        A2 = self.A2
        S = self.S
        y0 = ycenter
        x0 = xcenter

        values[0] = 0.
        values[1] = 0.
        values[2] = self.A1  + self.A2*exp(-((x[1]-y0)**2.+(x[0]-x0)**2.)/(2*S*S))
        values[3] = 0.

    def value_shape(self):
      return (4,)

# Inflow boundary
class InflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] < xmin + bmarg

# No-slip boundary
class NoSlipBoundary(SubDomain):
    def inside(self, x, on_boundary):
        dx = x[0] - xcenter
        dy = x[1] - ycenter
        r = sqrt(dx*dx + dy*dy)
        return on_boundary and \
               (x[1] < ymin + bmarg or x[1] > ymax - bmarg or \
                r < radius + bmarg)

class CylinderBoundary(SubDomain):
    def inside(self, x, on_boundary):
        r = sqrt((x[0]-xcenter)**2 + (x[1]-ycenter)**2)
        return on_boundary and r < radius + bmarg

# Outflow boundary
class OutflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] > xmax - bmarg

# Problem definition
class Problem(ProblemBase):

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Load mesh
        self.Nx = options["Nx"]
        self.U = Expression(('4*Um*x[1]*(H - x[1])/(H*H)', '0.0'), Um=Um, H=ymax)

        rect = Rectangle(xmin, ymin, xmax, ymax)
        circ = Circle(xcenter, ycenter, radius)
        domain = rect - circ
        self.mesh = Mesh(domain, self.Nx)

        #rescale Reynolds number to the problem
        Ubar = 2.*self.U((0,ymax/2.))[0]/3.
        options['Re'] = Ubar*Diameter*options['Re']

        self.t0 = 0.
        self.T = options['T']
        self.k = options['dt']

    def initial_conditions(self, W):
        w0 = InitialConditions()
        w0 = project(w0,W)

        return w0

    def boundary_conditions(self, W, t):
        # Create inflow boundary condition
        g0 = self.U
        f0 = Expression('A1  + A2*exp(-(pow(x[1]-y0,2)+pow(x[0]-x0,2))/(2*S*S))', \
                        A1=rhoMin, A2=rhoMax, S=S, y0=ycenter, x0=xcenter, t=t)

        bc0 = DirichletBC(W.sub(0), g0, InflowBoundary())

        # Create no-slip boundary condition
        bc1 = DirichletBC(W.sub(0), Constant((0,0)), NoSlipBoundary())

        # Create outflow boundary condition for pressure
        bc2 = DirichletBC(W.sub(2), Constant(0), OutflowBoundary())

        # Density boundary conditions for cool effects
        bc3 = DirichletBC(W.sub(1), rhoMin, InflowBoundary())
        bc4 = DirichletBC(W.sub(1), rhoMin, NoSlipBoundary())
        bc5 = DirichletBC(W.sub(1), f0, CylinderBoundary())

        # Collect boundary conditions
        bcs = [bc0, bc1, bc2, bc3, bc4, bc5]

        return bcs

    def F1(self, t):
        #forcing function for the momentum equation
        return Constant((0,0))

    def F2(self, t):
        #mass source for the continuity equation
        return Constant(0)

    def __str__(self):
        return "DensityCylinder"
