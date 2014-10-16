__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__ = "GNU GPL version 3 or any later version"

from AFES import *
from AFES import Solver as SolverBase


class Solver(SolverBase):

    '''
        Solver class for solving the Advection-Diffusion-Reaction equation.
    '''

    def __init__(self, options):

        try:
            self.epsilon = options['epsilon']
        except:
            self.epsilon = 1.
        try:
            self.a = options['a']
        except:
            self.a = 1
        try:
            self.beta = options['beta']
        except:
            self.beta = 1

        SolverBase.__init__(self, options)

    # strong residual for cG(1)cG(1)
    def strong_residual(self, u):
        a = self.a  # reaction coefficient
        beta = self.beta  # velocity

        R = a * u + dot(beta, grad(u))

        return R

    # weak residual for cG(1)cG(1)
    def weak_residual(self, problem, k, V, u, U, U_, v, ei_mode=False):
        h = CellSize(W.mesh())  # mesh size
        d = self.stabilization_parameters(U_, k, h)  # stabilization parameters

        a = self.a  # reaction coefficient
        epsilon = self.epsilon  # diffusion coefficient
        beta = self.beta  # velocity

        t0 = problem.t0

        t = t0 + k
        # forcing and mass source/sink
        F = problem.F(t)

        # least squares stabilization
        if ei_mode:
            d = 0

        # weak form of the equations
        r = ((1. / k) * inner(U - U_, v)
             + a * inner(u, v)
             + epsilon * inner(grad(u), grad(v))
             + inner(dot(beta, grad(u)), v)) * dx

        # forcing function
        r -= inner(F, v) * dx

        R = self.strong_residual(u)
        Rv = self.strong_residual(v)
        r += d * inner(R - F, Rv) * dx

        return r

    def stabilization_parameters(self, u, k, h):
        K = 1.
        if(h > self.epsilon):
            d = K * (k ** (-2) + inner(u, u) * h ** (-4)) ** (-0.5)
        else:
            d = K * (k ** (-2) + inner(u, u) * h ** (-2)) ** (-0.5)

        return d

    def __str__(self):
        return 'ADR'
