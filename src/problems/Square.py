__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__ = "GNU GPL version 3 or any later version"

from AFES import *
from AFES import Problem as ProblemBase


class Problem(ProblemBase):

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Create mesh
        Nx = options["Nx"]
        Ny = options["Ny"]

        try:
            x0 = float(options["x0"])
            x1 = float(options["x1"])
            y0 = float(options["y0"])
            y1 = float(options["y1"])
            self.mesh = RectangleMesh(x0, y0, x1, y1, Nx, Ny)
        except:
            self.mesh = UnitSquareMesh(Nx, Ny)

        self.t0 = 0.
        self.T = options['T']
        self.k = options['k']

        try:
            self.kappa = Constant(options['kappa'])
        except:
            self.kappa = Constant(1E-2)  # heat Coefficient
        try:
            self.rho = Constant(options['rho'])
        except:
            self.rho = Constant(1.)  # density

        try:
            self.c = Constant(options['c'])
        except:
            self.c = Constant(1.)  # heat Coefficient

    def initial_conditions(self, V):
        u0 = project(Expression('sin(pi*x[0])*sin(pi*x[1])'), V)

        return u0

    def boundary_conditions(self, V, t):
        # Create no-slip boundary condition for velocity
        bcs = DirichletBC(W, Constant(0.0), 'on_boundary')

        return bcs

    def F(self, t):
        # forcing function for the momentum equation
        return Constant(0)

    def __str__(self):
        return 'Square'
