__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__ = "GNU GPL version 3 or any later version"
#
#   adapted from channel.py in nsbench originally developed
#   by Kent-Andre Mardal <kent-and@simula.no>
#

from AFES import *
from AFES import Problem as ProblemBase

x0 = 0
x1 = 6
y0 = 0
y1 = 2


class InitialConditions(Expression):

    def eval(self, values, x):
        values[0] = 0.
        values[1] = 0.
        values[2] = 0.

    def value_shape(self):
        return (3,)


class InflowBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], x0)


class OutflowBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], x1)


class NoslipBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary and (near(x[1], y0) or near(x[1], y1))


class Problem(ProblemBase):

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Create mesh
        Nx = options["Nx"]
        Ny = options["Ny"]
        self.mesh = RectangleMesh(x0, y0, x1, y1, Nx, Ny)

        self.t0 = 0.
        self.T = options['T']
        self.k = options['k']

        try:
            self.nu = options['nu']
        except:
            self.nu = 1E-3

        self.Re = (y1 - y0)/self.nu

    def initial_conditions(self, W):
        w0 = InitialConditions()
        w0 = project(w0, W)

        return w0

    def boundary_conditions(self, W, t):

        noslip = DirichletBC(W.sub(0), Constant((0, 0)), NoslipBoundary())
        inflow = DirichletBC(W.sub(0), Constant((1, 0)), InflowBoundary())
        outflow = DirichletBC(W.sub(1), Constant(0), OutflowBoundary())

        bcs = [noslip, inflow, outflow]

        return bcs

    def F1(self, t):
        # forcing function for the momentum equation
        return Constant((0, 0))

    def F2(self, t):
        # mass source for the continuity equation
        return Constant(0)

    def functional(self, mesh, w):
        u, p = (as_vector((w[0], w[1])), w[2])

        M = u[0] * dx

        return M

    def __str__(self):
        return 'Channel'
