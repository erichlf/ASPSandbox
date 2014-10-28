j__author__ = "Erich L Foster <erichlf@gmail.com>"
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
import sys

geps = 1.0e-8

d = 1.
eta = 0.1
g = 9.8

x0 = -2 * d
x1 = 2 * d
y0 = -d
y1 = d

rhoMin = 1.
rhoMax = 1000
At = (rhoMax - rhoMin) / (rhoMin + rhoMax)  # Atwood number


class InitialConditions(Expression):

    def eval(self, values, x):
        values[0] = 0.
        values[1] = 0.
        if x[1] <= -0.7 * d:
            values[2] = rhoMax
        else:
            values[2] = rhoMin

        if x[0] <= -1.0 * d and x[1] <= -0.2 * d and x[0] >= -1.8 * d \
                and x[1] >= -1.0 * d:
            values[2] = rhoMax

        values[3] = 0.

    def value_shape(self):
        return (4,)


class InnerBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary and (x[0] <= 2 * d - geps
                                and x[0] >= -2 * d + geps
                                and x[1] <= d - geps and x[1] >= - d + geps)


class NoslipBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary


class PressureBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary and (fabs(x[0]) >= 2 * d)


class PsiMarker(Expression):

    def eval(self, values, x):
        ib = InnerBoundary()

        if(ib.inside(x, True)):
            values[0] = 1.0
        else:
            values[0] = 0.0


class Problem(ProblemBase):

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        try:
            self.mesh = Mesh(options['initialMesh'])
        except:
            print "ERROR: This problem requires a mesh. Please provide a mesh."
            sys.exit(1)

        self.psimarker = PsiMarker()

        self.t0 = 0.
        self.T = options['T']
        C_CFL = 5.
        self.Ubar = 1.
        self.k = C_CFL * self.mesh.hmin() / self.Ubar  # time step

    def initial_conditions(self, W):
        w0 = InitialConditions()
        ww0 = Function(W)
        w = TestFunction(W)
        F = inner(ww0, w) * dx + 1e-3 * inner(grad(ww0), grad(w)) * dx \
            - inner(w0, w) * dx
        solve(F == 0, ww0)

        return ww0

    def boundary_conditions(self, W, t):
        # Create no-slip boundary condition for velocity
        bc = DirichletBC(W.sub(0), Constant((0.0, 0.0)), NoslipBoundary())

        return bc

    def F(self, t):
        # forcing function for the momentum equation

        return Expression(('0.0', '-g'), g=g, t=t)

    def functional(self, W, w):
        (u, rho, p) = (as_vector((w[0], w[1])), w[2], w[3])
        # beta = Expression("1.0")
        # beta = Expression(('50.0*exp(-50.0*(pow(x[0] - 1.0, 2) '
        #                   + ' + pow(x[1] - -0.5, 2)))'))
        n = FacetNormal(W.mesh())
        # M = rho*beta*dx  # Mean of the x-velocity in the whole domain
        M = problem.psimarker*p*n[0]*ds  # Drag (only pressure)

        return M

    def __str__(self):
        return 'Dam'
