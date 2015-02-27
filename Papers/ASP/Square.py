__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__ = "GNU GPL version 3 or any later version"

from ASP import *
from ASP import Problem as ProblemBase


class Problem(ProblemBase):

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Create mesh
        Nx = options["Nx"]
        Ny = options["Ny"]

        self.mesh = UnitSquareMesh(Nx, Ny)

        self.t0 = 0.  # initial time
        self.T = options['T']  # final time
        self.k = options['k']  # time step

        self.kappa = Constant(1E-2)  # heat Coefficient

    def initial_conditions(self, V):
        u0 = project(Expression('sin(pi*x[0])*sin(pi*x[1])'), V)

        return u0

    def boundary_conditions(self, V, t):
        # Create no-slip boundary condition for velocity
        bcs = DirichletBC(V, Constant(0.0), 'on_boundary')

        return bcs

    def F(self, t):
        # forcing function for the momentum equation
        return Constant(0)

    def __str__(self):
        return 'Square'
