__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__ = "GNU GPL version 3 or any later version"

from AFES import *
from AFES import Problem as ProblemBase

# outer square dimensions
Xmin = 0.
Xmax = 1.
Ymin = 0.
Ymax = 1.
# inner square dimensions
xmin = 0.48
xmax = 0.52
ymin = 0.48
ymax = 0.52

kappa = Constant(1E-5)
a = Constant(1.)
beta = Constant((-1, -0.61))
rho = 1.
c = 1.

TR = 0.5
TA = 0.5
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
        self.Nx = options["Nx"]

        outerRect = Rectangle(Xmin, Ymin, Xmax, Ymax)
        innerRect = Rectangle(xmin, ymin, xmax, ymax)
        domain = outerRect - innerRect
        self.mesh = Mesh(domain, self.Nx)

        self.t0 = 0.
        self.T = options['T']
        self.k = options['k']

        try:  # set up heat coefficient
            self.kappa = Expression('kappa', kappa=options['kappa'])
        except:
            self.kappa = kappa
        try:  # velocity for adr
            self.beta = Expression('beta', beta=options['beta'])
        except:
            self.beta = beta
        try:  # reaction coefficient for adr
            self.a = Expression('a', a=options['a'])
        except:
            self.a = a

        self.rho = rho  # density
        self.c = c

    def initial_conditions(self, V):
        u0 = InitialConditions()
        u0 = project(u0, V)

        return u0

    def boundary_conditions(self, V, t):
        g = Expression('TR - TA*cos(omega*t)', TR=TR, TA=TA, omega=omega, t=t)
        bc0 = DirichletBC(V, Constant(0.0), OuterBoundary())
        bc1 = DirichletBC(V, g, InnerBoundary())

        return [bc0, bc1]

    def F(self, t):
        return Constant(0)

    def update(self, V, t):

        return self.boundary_conditions(V, t)

    def functional(self, V, u):

        psi = project(Expression('exp(-20 * (pow(x[0] - 0.25, 2) '
                                 + ' + pow(x[1] - 0.25, 2)))'), V)
        M = psi * dx

        return M

    def __str__(self):
        return 'Hole'
