__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"
#
#   adapted from channel.py in nsbench originally developed 
#   by Kent-Andre Mardal <kent-and@simula.no>
#

from problembase import *
from numpy import array

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
        N = options["N"]
        self.mesh = RectangleMesh(-1,-1,1,1,N, N)

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
        eta0 = Expression('A*exp(-(x[0]*x[0]+x[1]*x[1])/(2*S*S))', A=1.0, S=5E-2)

        return u0, eta0

    def boundary_conditions(self, V, Q, t):
        # Create no-slip boundary condition for velocity
        bcs = DirichletBC(V, Constant((0.0, 0.0)), NoslipBoundary())

        return bcs

    def __str__(self):
        return 'Drop'
