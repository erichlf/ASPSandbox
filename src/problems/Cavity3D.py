__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"
#
#   adapted from channel.py in nsbench originally developed 
#   by Kent-Andre Mardal <kent-and@simula.no>
#

from problembase import *
from numpy import array

# Inflow boundary
class Lid(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[2], 1.0)

# No-slip boundary
class NoslipBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (near(x[0], 0.0) or near(x[0], 1.0) \
                or near(x[1], 0.0) or near(x[1], 1.0) \
                or near(x[2], 0.0))

class Pressure(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (near(x[0], 0) or near(x[1], 0) or near(x[2], 0))

# Problem definition
class Problem(ProblemBase):
#   2D channel flow.

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Create mesh
        Nx = options["Nx"]
        Ny = options["Ny"]
        Nz = options["Nz"]
        self.mesh = UnitCubeMesh(Nx, Ny, Nz)

        self.t0 = 0.
        self.T = options['T']
        self.k = options['dt']

    def initial_conditions(self, V, Q):
        u0 = Constant((0, 0, 0))
        p0 = Constant(0)

        return u0, p0

    def boundary_conditions(self, W, t):
        # Create no-slip boundary condition for velocity
        noslip = DirichletBC(W.sub(0), Constant((0, 0, 0)), NoslipBoundary())

        # Create boundary conditions for pressure
        lid = DirichletBC(W.sub(0), Constant((1.0, 0, 0)), Lid())

        #pressure boundary
        pressure = DirichletBC(W.sub(1), Constant(0), Pressure())

        bcs = [noslip, lid, pressure]

        return bcs

    def F1(self, t):
        #forcing function for the momentum equation
        return Constant((0, 0, 0))

    def F2(self, t):
        #mass source for the continuity equation
        return Constant(0)

    def __str__(self):
        return 'Cavity3D'
