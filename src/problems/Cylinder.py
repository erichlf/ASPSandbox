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
        N = options["N"]
        self.mesh = Mesh('data/cylinder_0.xml.gz')

        self.nu = options["nu"] #viscosity
        self.rho = options["rho"] #density
        self.f0 = options["f0"] #reference Coriolis parameter
        self.beta = options["beta"] #beta plane parameter
        self.g = options["g"] #gravity
        self.h = options["h"] #fluid depth
        self.T = options["T"] #final time
        self.dt = options["dt"] #time-step
        self.theta = options["theta"] #theta for the theta time stepping method
        self.solver = options["linear_solver"] #what linear solver to use
        self.Pu = options["velocity_order"] #order of velocity element
        self.Pp = options["height_order"] #order of height/pressure element
        self.stabilize = options["stabilize"] #should we stabilize the solver 

    def initial_conditions(self, V, Q):

        u0 = Constant((0, 0))
        p0 = Constant(0)

        return u0, p0

    def boundary_conditions(self, V, Q, t):

        # Create inflow boundary condition
        self.g0 = Expression(('4*Um*(x[1]*(ymax-x[1]))*sin(pi*t/8.0)/(ymax*ymax)', '0.0'),
                             Um=1.5, ymax=ymax, t=t)
        self.b0 = InflowBoundary()
        bc0 = DirichletBC(V, self.g0, self.b0)

        # Create no-slip boundary condition
        self.b1 = NoslipBoundary()
        self.g1 = Constant((0, 0))
        bc1     = DirichletBC(V, self.g1, self.b1)

        # Create outflow boundary condition for pressure
        self.b2 = OutflowBoundary()
        self.g2 = Constant(0)
        bc2     = DirichletBC(Q, self.g2, self.b2)

        # Collect boundary conditions
        bcs = [bc0, bc1, bc2]

        return bcs

    def update(self, t, u, p):
        self.g0.t = t

    def F(self, t):
        return Constant((0,0))

    def __str__(self):
        return "Cylinder"
