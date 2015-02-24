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

    def function_space(self, mesh):  # define functions spaces
        try:
            Pu = options['velocity_order']
        except:
            Pu = 1
        try:
            Pp = options['pressure_order']
        except:
            Pp = 1

        if Pu == Pp and Pu == 1:
            self.stabilize = True
        else:
            self.stabilize = False

        V = VectorFunctionSpace(mesh, 'CG', Pu)
        Q = FunctionSpace(mesh, 'CG', Pp)
        W = MixedFunctionSpace([V, Q])

        return W

    # strong residual for cG(1)cG(1)
    def strong_residual(self, W, w, w2):
        if W.mesh().topology().dim() == 2:
            (u, p) = (as_vector((w[0], w[1])), w[2])
            u2 = as_vector((w2[0], w2[1]))
        else:
            (u, p) = (as_vector((w[0], w[1], w[2])), w[3])
            u2 = as_vector((w2[0], w2[1], w2[2]))

        R1 = grad(u2) * u + grad(p)
        R2 = div(u)

        return R1, R2

    # weak residual for cG(1)cG(1)
    def weak_residual(self, problem, k, W, w_theta, w, w_, wt, ei_mode=False):

        if W.mesh().topology().dim() == 2:
            (u, p) = (as_vector((w_theta[0], w_theta[1])), w_theta[2])
            U = as_vector((w[0], w[1]))
            U_ = as_vector((w_[0], w_[1]))
            (v, q) = (as_vector((wt[0], wt[1])), wt[2])
        else:
            (u, p) = (as_vector((w_theta[0], w_theta[1], w_theta[2])),
                      w_theta[3])
            U = as_vector((w[0], w[1], w[2]))
            U_ = as_vector((w_[0], w_[1], w_[2]))
            (v, q) = (as_vector((wt[0], wt[1]), wt[2]), wt[3])

        nu = problem.nu

        h = CellSize(W.mesh())
        d1 = conditional(le(h, nu), h**2, h)  # stabilization parameter

        if ei_mode or not self.stabilize:  # turn off stabilization in ei_mode
            d1 = Constant(0)

        f = problem.F1(problem.t0)  # forcing

        # Weak form
        r = (1. / k * inner(U - U_, v)
             + nu*inner(grad(u), grad(v))
             + inner(grad(p) + grad(u)*u, v)
             + div(u)*q)*dx
        r -= inner(f, v)*dx

        # GLS stabilization
        R1, R2 = self.strong_residual(W, w_theta, w_theta)
        Rv1, Rv2 = self.strong_residual(W, wt, w_theta)
        r += d1*(inner(R1 - f, Rv1) + R2 * Rv2) * dx

        return r

    def stabilization_parameters(self, u, p, k, h):
        K1 = 1.
        K2 = 1.
        d1 = K1 * h  # (k**(-2) + inner(u,u)*h**(-2))**(-0.5)
        d2 = K2 * h

        return d1, d2

    def suffix(self, problem):
        '''
            Obtains the run specific data for file naming, e.g. Nx, k, etc.
        '''

        try:
            s = 'Re' + str(problem.Re)
        except:
            s = 'nu' + str(problem.nu)

        s += 'T' + str(problem.T)
        if problem.Nx is not None:
            s += 'Nx' + str(problem.Nx)
        if problem.Ny is not None:
            s += 'Ny' + str(problem.Ny)
        if problem.mesh.topology().dim() > 2 and problem.Nz is not None:
            s += 'Nz' + str(problem.Nz)

        s += 'K' + str(problem.k)

        return s

    def __str__(self):
        return 'NSE'
