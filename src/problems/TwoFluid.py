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

from problembase import *

d = 1.
eta = 0.1
g = 9.8

x0 = -d / 2.
x1 = d / 2.
y0 = -2 * d
y1 = 2 * d

At = 0.5  # Atwood number
rhoMin = 1.
rhoMax = rhoMin * (1. + At) / (1. - At)


class InitialConditions(Expression):

    def eval(self, values, x):
        values[0] = 0.
        values[1] = 0.
        values[2] = 0.5 * (rhoMin + rhoMax) + 0.5 * (rhoMax - rhoMin) * \
            tanh((x[1] + eta * cos(2 * pi * x[0] / d)) / (0.01 * d))
        values[3] = 0.

    def value_shape(self):
        return (4,)


class NoslipBoundary(SubDomain):
    # No-slip boundary

    def inside(self, x, on_boundary):
        return on_boundary


class Problem(ProblemBase):
    # Problem definition

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Create mesh
        self.Nx = options['Nx']
        self.Ny = options['Ny']
        self.mesh = RectangleMesh(x0, y0, x1, y1, self.Nx, self.Ny)
        self.Fr = options['Fr']

        self.t0 = 0.
        self.T = options['T']
        self.k = options['dt']

    # get the initial condition and project it
    def initial_conditions(self, W):
        w0 = InitialConditions()
        w0 = project(w0, W)

        return w0

    def boundary_conditions(self, W, t):
        # Create no-slip boundary condition for velocity
        bcs = DirichletBC(W.sub(0), Constant((0.0, 0.0)), NoslipBoundary())

        return bcs

    def F1(self, t):
        # forcing function for the momentum equation
        return Expression(('0.', '1./g'), g=g, t=t)

    def F2(self, t):
        # mass source for the continuity equation
        return Expression('0.0', t=t)

    def functional(self, mesh, w):

        (u, rho, p) = (as_vector((w[0], w[1])), w[2], w[3])

        M = rho * dx  # Mean of the x-velocity in the whole domain

        return M

    def __str__(self):
        return 'TwoFluid'
