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
        self.nu = 1. / options['Re']

    # strong residual for cG(1)cG(1)
    def strong_residual(self, W, w, w2):
        if W.mesh().topology().dim() == 2:
            (u, p) = (as_vector((w[0], w[1])), w[2])
            (U, P) = (as_vector((w2[0], w2[1])), w2[2])
        else:
            (u, p) = (as_vector((w[0], w[1], w[2])), w[3])
            (U, P) = (as_vector((w2[0], w2[1], w2[2])), w2[3])

        R1 = grad(U) * u + grad(p)
        R2 = div(u)

        return R1, R2

    # weak residual for cG(1)cG(1)
    def weak_residual(self, problem, k, W, w, ww, w_, wt, ei_mode=False):
        if W.mesh().topology().dim() == 2:
            (u, p) = (as_vector((w[0], w[1])), w[2])
            (U, P) = (as_vector((ww[0], ww[1])), ww[2])
            (U_, P_) = (as_vector((w_[0], w_[1])), w_[2])
            (v, q) = (as_vector((wt[0], wt[1])), wt[2])
        else:
            (u, p) = (as_vector((w[0], w[1], w[2])), w[3])
            (U, P) = (as_vector((ww[0], ww[1], ww[2])), ww[3])
            (U_, P_) = (as_vector((w_[0], w_[1], w_[2])), w_[3])
            (v, q) = (as_vector((wt[0], wt[1], wt[2])), wt[3])

        h = CellSize(W.mesh())  # mesh size
        d1, d2 = self.stabilization_parameters(
            U_, P_, k, h)  # stabilization parameters

        nu = self.nu  # Reynolds Number

        t0 = problem.t0

        t = t0 + k
        # forcing and mass source/sink
        F = problem.F1(t)

        # least squares stabilization
        if(ei_mode):
            d1 = 0
            d2 = 0

        # weak form of the equations
        r = (1. / k) * inner(U - U_, v) * dx \
            + inner(grad(p) + grad(u) * u, v) * dx
        r += nu * inner(grad(u), grad(v)) * dx
        r += div(u) * q * dx

        # forcing function
        r -= inner(F, v) * dx

        R1, R2 = self.strong_residual(W, w, w)
        Rv1, Rv2 = self.strong_residual(W, wt, w)
        r += (d1 * inner(R1 - F, Rv1) + d2 * R2 * Rv2) * dx

        return r

    def stabilization_parameters(self, u, p, k, h):
        K1 = 1.
        K2 = 1.
        d1 = K1 * h  # (k**(-2) + inner(u,u)*h**(-2))**(-0.5)
        d2 = K2 * h

        return d1, d2

    def __str__(self):
        return 'NSE'
