__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__ = "GNU GPL version 3 or any later version"

from AFES import *
from AFES import Solver as SolverBase


class Solver(SolverBase):

    '''
        Solver class for solving both 2D and 3D Incompressible Navier-Stokes
        equation.
    '''

    def __init__(self, options):
        SolverBase.__init__(self, options)

    # strong residual for cG(1)cG(1)
    def strong_residual(self, W, w, w2):
        (u, p) = (as_vector((w[0], w[1])), w[2])

        rho = self.rho  # density
        c = self.c
        Ud = self.Ud

        R1 = -rho * grad(u) * Ud + grad(p)
        R2 = -1. / (rho * c * c) * dot(Ud, grad(p)) + div(u)

        return R1, R2

    # weak residual for cG(1)cG(1)
    def weak_residual(self, problem, k, W, w_theta, w, w_, wt, ei_mode=False):
        (u, p) = (as_vector((w_theta[0], w_theta[1])), w_theta[2])
        (U, P) = (as_vector((w[0], w[1])), w[2])
        (U_, P_) = (as_vector((w_[0], w_[1])), w_[2])
        (v, q) = (as_vector((wt[0], wt[1])), wt[2])

        h = CellSize(W.mesh())  # mesh size

        rho = problem.rho  # density
        c = problem.c  # speed of sound
        self.rho = rho
        self.c = c

        Ud = problem.Ud
        self.Ud = Ud

        t0 = problem.t0

        t = t0 + k
        # forcing and mass source/sink
        F = problem.F(t)
        Q = problem.Q(t)

        d1, d2 = self.stabilization_parameters(h)

        # weak residual for cG(1)cG(1)
        r = (rho / k) * inner(U - U_, v) * dx \
            + inner(grad(p) + rho * grad(u) * Ud, v) * dx
        r += 1. / (rho * c * c) * (P - P_) * q * dx \
            - (1. / (rho * c * c) * dot(Ud, grad(p)) - div(u)) * q * dx

        # forcing function
        r -= (inner(F, v) + Q * q) * dx

        # least squares stabilization
        R1, R2 = self.strong_residual(W, w_theta, w_theta)
        Rv1, Rv2 = self.strong_residual(W, wt, w_theta)
        r += (d1 * inner(R1 - F, Rv1) + d2 * (R2 - Q) * Rv2) * dx

        return r

    def stabilization_parameters(self, h):
        K1 = 1.
        K2 = 1.
        d1 = K1 * h
        d2 = K2 * h

        return d1, d2

    def suffix(self, problem):
        import numpy as np
        '''
            Obtains the run specific data for file naming, e.g. Nx, k, etc.
        '''
        s = ''

        Ud = np.zeros(2)
        problem.Ud.eval(Ud, np.zeros(3))
        Ud = np.dot(Ud, Ud)**0.5

        # Return file suffix for output files
        s = 'rho%.3Gc%.3GUd%.3G' % (problem.rho, problem.c, Ud)

        s += 'T%G' % (problem.T)
        if problem.Nx is not None:
            s += 'Nx' + str(problem.Nx)
        if problem.Ny is not None:
            s += 'Ny' + str(problem.Ny)
        if problem.mesh.topology().dim() > 2 and problem.Nz is not None:
            s += 'Nz' + str(problem.Nz)

        s += 'K' + str(problem.k)

        return s

    def __str__(self):
        return 'MixdWaveALE'
