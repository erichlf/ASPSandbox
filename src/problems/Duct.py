__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2014-10-14"
__license__ = "GNU GPL version 3 or any later version"
#
#   adapted from channel.py in nsbench originally developed
#   by Kent-Andre Mardal <kent-and@simula.no>
#

from ASP import *
from ASP import Problem as ProblemBase

# domain dimensions
x0 = 0.
x1 = 17.
y0 = 0.
y1 = 2.

rho = 1.225  # density of air in kg/m^3
c = 350  # speed of sound in m/s

n = 3  # nth resonance
f = (2 * n - 1) * c / (4 * (x1 - x0))  # nth resonance frequency
Ud = Expression(('0', '0'))  # domain velocity


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


class NoNormalBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary and (near(x[1], y0) or near(x[1], y1))


class Problem(ProblemBase):

    '''
        Here we define a long duct for the Mixed Wave Equation
    '''

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Create mesh
        Nx = options["Nx"]
        Ny = options["Ny"]
        self.mesh = RectangleMesh(x0, y0, x1, y1, Nx, Ny)

        self.t0 = 0.
        self.T = options['T']
        self.k = options['k']

        self.rho = rho  # density
        self.c = c  # speed of sound

        self.Ug = Expression(('sin(2*pi*f*t)', '0'), f=f, t=self.t0)  # BC
        self.Ud = Ud  # domain velocity

    def initial_conditions(self, W):
        w0 = InitialConditions()

        return project(w0, W)

    def boundary_conditions(self, W, t):
        self.Ug.t = t

        # Create no-normal flow boundary condition for velocity
        noNormal = DirichletBC(W.sub(0).sub(1), Constant(0), NoNormalBoundary())

        inflow = DirichletBC(W.sub(0), self.Ug, InflowBoundary())

        # Create boundary conditions for pressure
        outflow = DirichletBC(W.sub(1), Constant(0), OutflowBoundary())

        bcs = [noNormal, inflow, outflow]

        return bcs

    def update(self, W, t):

        return self.boundary_conditions(W, t)

    def F(self, t):  # forcing function for the momentum equation
        return Constant((0, 0))

    def Q(self, t):  # mass source for the continuity equation
        return Constant(0)

    def __str__(self):
        return 'Duct'
