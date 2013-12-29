__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"
#
#   adapted from channel.py in nsbench originally developed 
#   by Kent-Andre Mardal <kent-and@simula.no>
#

from problembase import *
from numpy import array

# Inflow boundary
class InflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0],0.)

# Outflow boundary
class OutflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0],1.)

# No-slip boundary
class NoslipBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (near(x[1],0.) or near(x[1], 1.))

# Problem definition
class Problem(ProblemBase):
#   2D channel flow.

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Create mesh
        N = options["N"]
        self.mesh = UnitSquare(N, N)

    def initial_conditions(self, V, Q):
        u0 = Constant((0, 0))
        p0 = Constant(0)#self.pressure_bc(V, Q, 0)

        return u0, p0

    def boundary_conditions(self, V, Q, t):
        # Create no-slip boundary condition for velocity
        noslip = DirichletBC(V, (0, 0), NoslipBoundary())

        # Create boundary conditions for pressure
        inflow = DirichletBC(V, Constant((1,0)), InflowBoundary())
        outflow = DirichletBC(Q, Constant(0), OutflowBoundary())

        bcs = [noslip, inflow, outflow]

        return bcs

    def pressure_bc(self, V, Q, t):
        return Expression('p*(1.-x[0])',p=1.,t=t)

    def F1(self, t):
        #forcing function for the momentum equation
        return Constant((0,0))

    def F2(self, t):
        #mass source for the continuity equation
        return Constant(0)

    def __str__(self):
        return 'Channel'
