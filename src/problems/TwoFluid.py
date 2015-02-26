__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__ = "GNU GPL version 3 or any later version"

from AFES import *
from AFES import Problem as ProblemBase

d = 1.
eta = 0.01
g = 10

x0, x1 = -d / 2., d / 2.
y0, y1 = -2 * d, 2 * d


class InitialConditions(Expression):

    def __init__(self, At):
        self.rhoMin = 1.
        self.rhoMax = self.rhoMin * (1. + At) / (1. - At)

    def eval(self, values, x):
        rhoMin = self.rhoMin
        rhoMax = self.rhoMax

        values[0], values[1], values[3] = 0., 0., 0.
        values[2] = 0.5 * (rhoMin + rhoMax) + 0.5 * (rhoMax - rhoMin) * \
            tanh((x[1] - 1.5 * d + eta * cos(2 * pi * x[0] / d)) / (0.01 * d))

    def value_shape(self):
        return (4,)


class NoslipBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary


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
        self.Nx, self.Ny = options['Nx'], options['Ny']
        self.mesh = RectangleMesh(x0, y0, x1, y1, self.Nx, self.Ny)

        self.t0, self.T, = 0., options['T']
        try:
            self.CFL = options['CFL']
        except:
            self.CFL = 5.
        self.Ubar = g
        self.k = self.time_step(self.Ubar, self.mesh)  # mesh size

    def time_step(self, u, mesh):
        return self.CFL * mesh.hmin()/u

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

        F = (wb * wt + h * inner(grad(wb), grad(wt))) * dx

        solve(F == 0, wb, bc0)
        self.wb = wb

    def boundary_conditions(self, W, t):

        bc1 = DirichletBC(W.sub(0), Constant((0.0, 0.0)), NoslipBoundary())

        return [bc1]

    def F(self, t):
        # forcing function for the momentum equation
        return Expression(('0.', '-g'), g=g, t=t)

    def functional(self, W, w):

        u = as_vector((w[0], w[1]))

        return (u[1] + u[0]) * dx

    def __str__(self):
        return 'TwoFluid'
