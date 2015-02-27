__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"

from ASP import *
from ASP import Solver as SolverBase


class Solver(SolverBase):

    '''
        Solver class for solving the Heat equation.
    '''

    def __init__(self, options):
        SolverBase.__init__(self, options)

    def function_space(self, mesh):
        # Define function spaces
        V = FunctionSpace(mesh, 'CG', self.Pu)

        return V

    def weak_residual(self, problem, k, V, u, U, U_, v, ei_mode=False):
        t = problem.t0

        # problem parameters
        kappa = problem.kappa  # heat coefficient
        rho = problem.rho  # density
        c = problem.c

        f = problem.F(t)  # forcing and mass source/sink

        # weak form of the equations
        r = rho * c * (1. / k) * (U - U_) * v * dx
        if kappa.rank() == 0:
            r += kappa * inner(grad(u), grad(v)) * dx
        else:  # anisotropic case
            r += inner(dot(kappa, grad(u)), grad(v)) * dx
        r -= f * v * dx

        return r

    def Plot(self, problem, V, u):

        # Plot velocity and height and wave object
        if self.vizU is None:
            self.vizU = plot(u, title='Temperature', rescale=True, elevate=0.0)
        else:
            self.vizU.plot(u)

    def __str__(self):
        return 'Heat'
