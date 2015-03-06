__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__ = "GNU GPL version 3 or any later version"

from ASP import *
from ASP import Problem as ProblemBase
from mshr import *

# outer square dimensions
Xmin, Xmax, Ymin, Ymax = 0., 1., 0., 1.
# inner square dimensions
xmin, xmax, ymin, ymax = 0.48, 0.52, 0.48, 0.52

kappa = Constant(1E-1)
alpha = Constant(1.)
beta = Constant((-1, -0.61))
rho, c = 1., 1.

TR, TA = 10., 10.
omega = 2. * pi


class InitialConditions(Expression):

    def eval(self, values, x):
        values[0] = 0.


class OuterBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary and (near(x[0], Xmin) or near(x[0], Xmax)
                                or near(x[1], Ymin) or near(x[1], Ymax))


class InnerBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary and (near(x[0], xmin) or near(x[0], xmax)
                                or near(x[1], ymin) or near(x[1], ymax))


# Problem definition
class Problem(ProblemBase):

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Create mesh

        try:
            domain = options['initial_mesh']
            self.mesh = Mesh(domain)
        except:
            outerRect = Rectangle(Point(Xmin, Ymin), Point(Xmax, Ymax))
            innerRect = Rectangle(Point(xmin, ymin), Point(xmax, ymax))
            domain = outerRect - innerRect
            self.mesh = generate_mesh(domain, self.Nx)

        self.t0 = 0.
        self.T = options['T']
        self.Ubar = 1.0
        self.k = self.time_step(self.Ubar, self.mesh)

        try:  # set up heat coefficient
            self.kappa = Constant(options['kappa'])
        except:
            self.kappa = kappa
        try:  # velocity for adr
            self.beta = Constant(options['beta'])
        except:
            self.beta = beta
        try:  # reaction coefficient for adr
            self.alpha = Constant(options['alpha'])
        except:
            self.alpha = alpha

        self.rho = rho  # density
        self.c = c

    def initial_conditions(self, V):
        u0 = InitialConditions()
        u0 = project(u0, V)

        return u0

    def boundary_conditions(self, V, t):
        F = Expression('TA/2 + TA/2*sin(pi*t)', TA=TA, t=t)
        bc0 = DirichletBC(V, Constant(0.0), OuterBoundary())
        bc1 = DirichletBC(V, F, InnerBoundary())

        return [bc0, bc1]

    def F(self, t):
        return Constant(0)

    def update(self, V, t):

        return self.boundary_conditions(V, t)

    def functional(self, V, u):

        psi = Expression('exp(-20 * (pow(x[0] - 0.25, 2) '
                         + ' + pow(x[1] - 0.25, 2)))')
        M = psi * u * dx

        return M

    def time_step(self, Ubar, mesh):
        CFL = 100

        k = CFL*mesh.hmin()

        return k

    def __str__(self):
        return 'Hole'
