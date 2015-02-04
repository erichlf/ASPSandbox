__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__ = "GNU GPL version 3 or any later version"
#
#   adapted from channel.py in nsbench originally developed
#   by Kent-Andre Mardal <kent-and@simula.no>
#

'''
This is an extremely boring problem with no forcing and on a square.
This is basically a blank problem that we can adapt with optional inputs.
'''

from AFES import *
from AFES import Problem as ProblemBase

d = 1.
eta = 0.01
g = 10

x0 = -d / 2.
x1 = d / 2.
y0 = -2 * d
y1 = 2 * d


class InitialConditions(Expression):

    def __init__(self, At):
        self.rhoMin = 1.
        self.rhoMax = self.rhoMin * (1. + At) / (1. - At)

    def eval(self, values, x):
        rhoMin = self.rhoMin
        rhoMax = self.rhoMax

        values[0] = 0.
        values[1] = 0.
        values[2] = 0.5 * (rhoMin + rhoMax) + 0.5 * (rhoMax - rhoMin) * \
            tanh((x[1] - 1.5 * d + eta * cos(2 * pi * x[0] / d)) / (0.01 * d))
        values[3] = 0.

    def value_shape(self):
        return (4,)


class NoslipBoundary(SubDomain):
    # No-slip boundary

    def inside(self, x, on_boundary):
        return on_boundary


class BottomBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary and x[1] < y1 - DOLFIN_EPS


class TopBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary and x[1] > y0 + DOLFIN_EPS


class Problem(ProblemBase):
    # Problem definition

    def __init__(self, options):

        ProblemBase.__init__(self, options)

        try:  # viscosity
            self.nu = options['nu']
        except:
            self.nu = 1E-3

        try:  # Atwood number
            self.At = options['At']
        except:
            self.At = 0.5

        self.Re = d**1.5*g/self.nu

        # Create mesh
        self.Nx = options['Nx']
        self.Ny = options['Ny']
        self.mesh = RectangleMesh(x0, y0, x1, y1, self.Nx, self.Ny)

        self.t0 = 0.
        self.T = options['T']
        self.k = options['k']

    def initial_conditions(self, W):

        # artificial viscosity for stabilization
        self.artificial_viscosity(W)

        w0 = InitialConditions(self.At)
        w0 = project(w0, W)

        return w0

    def artificial_viscosity(self, W):
        V = FunctionSpace(W.mesh(), 'CG', 1)

        bc0 = DirichletBC(V, Constant(1.0), NoslipBoundary())

        h = CellSize(V.mesh())

        wb = Function(V)
        wt = TestFunction(V)

        F = wb*wt*dx + h*inner(grad(wb), grad(wt))*dx

        solve(F == 0, wb, bc0)
        self.wb = wb

    def boundary_conditions(self, W, t):
        rhoMin = 1.
        rhoMax = rhoMin * (1. + self.At) / (1. - self.At)

        # Create no-slip boundary condition for velocity
        bc1 = DirichletBC(W.sub(0), Constant((0.0, 0.0)), NoslipBoundary())
        # bc2 = DirichletBC(W.sub(1), Constant(rhoMin), BottomBoundary())
        # bc3 = DirichletBC(W.sub(1), Constant(rhoMax), TopBoundary())

        return [bc1]  # , bc2, bc3]

    def F(self, t):
        # forcing function for the momentum equation
        return Expression(('0.', '-g'), g=g, t=t)

    def functional(self, W, w):

        (u, rho, p) = (as_vector((w[0], w[1])), w[2], w[3])

        M = (u[1] + u[0]) * dx

        return M

    def __str__(self):
        return 'TwoFluid'
