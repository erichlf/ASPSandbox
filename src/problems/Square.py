__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"
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
        eta0 = Constant(0)

        return u0, eta0

    def boundary_conditions(self, V, Q, t):
        # Create no-slip boundary condition for velocity
        bcs = DirichletBC(V, Constant((0.0, 0.0)), NoslipBoundary())

        return bcs

    def F(self, t):
        return Constant((0,0))

    def __str__(self):
        return 'Drop'
