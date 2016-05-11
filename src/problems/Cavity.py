__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__ = "GNU GPL version 3 or any later version"
#
#   adapted from channel.py in nsbench originally developed
#   by Kent-Andre Mardal <kent-and@simula.no>
#

from ASP import *
from ASP import Problem as ProblemBase
from mshr import *


class InitialConditions2D(Expression):

    def eval(self, values, x):
        values[0] = 0.
        values[1] = 0.
        values[2] = 0.

    def value_shape(self):
        return (3,)


class InitialConditions3D(Expression):

    def eval(self, values, x):
        values[0] = 0.
        values[1] = 0.
        values[2] = 0.
        values[3] = 0.

    def value_shape(self):
        return (4,)


class Lid(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 1.0)


class NoslipBoundary(SubDomain):

    def __init__(self, dim):
        SubDomain.__init__(self)
        self.dim = dim

    def inside(self, x, on_boundary):
        return on_boundary and (near(x[0], 0) or near(x[0], 1)
                                or near(x[1], 0.0)
                                or (self.dim == 3 and (near(x[1], 1)
                                                       or near(x[2], 0))))


class Pressure(SubDomain):

    def __init__(self, dim):
        SubDomain.__init__(self)
        self.dim = dim

    def inside(self, x, on_boundary):
        return on_boundary and (near(x[0], 0) or near(x[1], 0) or
                                (self.dim == 3 and near(x[2], 0)))


class Problem(ProblemBase):

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Create mesh
        self.dim = options['dim']
        self.Nx = options["Nx"]
        if self.dim == 3:
            self.uNoSlip = Expression(('0', '0', '0'))
            self.uLid = Expression(('1', '0', '0'))
            domain = Box(Point(0, 0, 0), Point(1, 1, 1))
        else:
            self.uNoSlip = Expression(('0', '0'))
            self.uLid = Expression(('1', '0'))
            domain = Rectangle(Point(0, 0), Point(1, 1))

        self.mesh = generate_mesh(domain, self.Nx)

        self.t0 = 0.
        self.T = options['T']
        self.k = options['k']

        try:
            self.nu = options['nu']
        except:
            self.nu = 1E-3

        self.Re = 1. / self.nu

    def initial_conditions(self, W, annotate=False):
        if self.dim == 2:
            w0 = InitialConditions2D()
        else:
            w0 = InitialConditions3D()

        w0 = project(w0, W)

        return w0

    def boundary_conditions(self, W, t):

        noslip = DirichletBC(W.sub(0), self.uNoSlip, NoslipBoundary(self.dim))
        lid = DirichletBC(W.sub(0), self.uLid, Lid())
        pressure = DirichletBC(W.sub(1), Constant(0), Pressure(self.dim))

        bcs = [noslip, lid, pressure]

        return bcs

    def F1(self, t):  # forcing function for the momentum equation
        if self.dim == 2:
            f = Constant((0, 0))
        else:
            f = Constant((0, 0, 0))

        return f

    def F2(self, t):  # mass source for the continuity equation
        return Constant(0)

    def functional(self, mesh, w):
        if self.dim == 2:
            u, p = (as_vector((w[0], w[1])), w[2])
        else:
            u, p = (as_vector((w[0], w[1], w[2])), w[3])

        M = u[0] * dx

        return M

    def __str__(self):
        return 'Cavity'
