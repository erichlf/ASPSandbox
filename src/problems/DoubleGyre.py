__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2014-07-17"
#
#   adapted from channel.py in nsbench originally developed 
#   by Kent-Andre Mardal <kent-and@simula.no>
#

from problembase import *
from numpy import array

x0 = 0
x1 = 1
y0 = -1
y1 = 1

# No-slip boundary
class NoslipBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary 

# Problem definition
class Problem(ProblemBase):
#  Double gyre forcing 

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Create mesh
        Nx = options['Nx']
        Ny = options['Ny']
        self.mesh = RectangleMesh(x0, y0, x1, y1, Nx, Ny)

        self.t0 = 0.
        self.T = options['T']
        self.k = options['dt']

    def initial_conditions(self, Q, Psi):
        q0 = Constant(0.)
        psi0 = Constant(0.)

        return q0, psi0

    def boundary_conditions(self, Q, Psi, t):
        # Create no-slip boundary condition for velocity
        noslipQ = DirichletBC(Q, Constant(0.), NoslipBoundary())
        noslipPsi = DirichletBC(Psi, Constant(0.), NoslipBoundary())

        bcs = [noslipQ, noslipPsi]

        return bcs

    def F(self, t):
        return Expression('sin(pi*x[1])', t=t) #Forcing function 

    def __str__(self):
        return 'DoubleGyre'
