__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2014-07-17"
#
#   adapted from channel.py in nsbench originally developed
#   by Kent-Andre Mardal <kent-and@simula.no>
#

from AFES import *
from AFES import Problem as ProblemBase

x0 = 0
x1 = 1
y0 = -1
y1 = 1


class InitialConditions(Expression):

    def eval(self, values, x):
        values[0] = 0.
        values[1] = 0.

    def value_shape(self):
        return (2,)


class NoslipBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary


class Problem(ProblemBase):

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Create mesh
        Nx = options['Nx']
        Ny = options['Ny']
        self.mesh = RectangleMesh(x0, y0, x1, y1, Nx, Ny)

        self.t0 = 0.
        self.T = options['T']
        self.k = options['k']

    def initial_conditions(self, W):
        w0 = InitialConditions()
        w0 = project(w0, W)

        return w0

    def boundary_conditions(self, W, t):
        # Create no-slip boundary condition for velocity
        noslipQ = DirichletBC(W.sub(0), Constant(0.), NoslipBoundary())
        noslipPsi = DirichletBC(W.sub(1), Constant(0.), NoslipBoundary())

        bcs = [noslipQ, noslipPsi]

        return bcs

    def F(self, t):
        return Expression('sin(pi*x[1])', t=t)  # Forcing function

    def functional(self, W, w):

        (q, psi) = (w[0], w[1])

        M = q * dx  # Mean of the vorticity in the whole domain

        return M

    def __str__(self):
        return 'DoubleGyre'
