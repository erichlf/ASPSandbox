__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"
#
#   adapted from channel.py in nsbench originally developed 
#   by Kent-Andre Mardal <kent-and@simula.no>
#

from problembase import *
from numpy import array

# Inflow boundary
class InflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS

# Outflow boundary
class OutflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > 1 - DOLFIN_EPS

# No-slip boundary
class NoslipBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS

# Problem definition
class Problem(ProblemBase):
#   2D channel flow.

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Create mesh
        N = options["N"]
        self.mesh = UnitSquare(N, N)

        self.F = Constant((0, 0)) #Forcing function 
        self.nu = options["nu"] #viscosity
        self.f0 = options["f0"] #reference Coriolis parameter
        self.beta = options["beta"] #beta plane parameter
        self.g = options["g"] #gravity
        self.h = options["h"] #fluid depth
        self.T = options["T"] #final time
        self.dt = options["dt"] #time-step

    def initial_conditions(self, V, Q):
        u0 = Constant((0, 0))
        p0 = Expression('1 - x[0]')

        return u0, p0

    def boundary_conditions(self, V, Q, t):
        # Create no-slip boundary condition for velocity
        bv = DirichletBC(V, Constant((0.0, 0.0)), NoslipBoundary())

        # Create boundary conditions for pressure
        bp0 = DirichletBC(Q, self.pressure_bc(Q), InflowBoundary())
        bp1 = DirichletBC(Q, self.pressure_bc(Q),  OutflowBoundary())

        bcs = [bv, bp0, bp1]

        return bcs

    def pressure_bc(self, Q):
        element = FiniteElement('CG', triangle, 1)
        return Expression('1 - x[0]', element=element)

    def __str__(self):
        return 'Channel'
