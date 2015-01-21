from AFES import *
from AFES import Solver as SolverBase


class Solver(SolverBase):
    def __init__(self, options):
        SolverBase.__init__(self, options)

        self.steady_state = True
        self.nu = options['nu']

    def strong_residual(self, w, w2):
        (u, p) = (as_vector((w[0], w[1])), w[2])
        (u2, p2) = (as_vector((w2[0], w2[1])), w2[2])

        R1 = grad(p) + grad(u2)*u
        R2 = div(u)

        return R1, R2

    def weak_residual(self, problem, W, w, wt, ei_mode=False):

        (u, p) = (as_vector((w[0], w[1])), w[2])
        (v, q) = (as_vector((wt[0], wt[1])), wt[2])

        nu = self.nu

        h = CellSize(W.mesh())
        d1 = conditional(le(h, nu), h**2, h)  # stabilization parameter

        if ei_mode:  # turn off stabilization in ei_mode
            d1 = Constant(0)

        f = problem.F()  # forcing

        # Weak form
        r = (nu*inner(grad(u), grad(v)) + inner(grad(p) + grad(u)*u, v)
             + div(u)*q)*dx
        r -= inner(f, v)*dx

        # GLS stabilization
        R1, R2 = self.strong_residual(w, w)
        Rv1, Rv2 = self.strong_residual(wt, w)
        r += d1*(inner(R1 - f, Rv1) + R2 * Rv2) * dx

        return r

    def suffix(self, problem):
        '''
            Obtains the run specific data for file naming, e.g. Nx, etc.
        '''

        try:
            s = 'Re' + str(problem.Re)
        except:
            s = 'nu' + str(problem.nu)

        if problem.Nx is not None:
            s += 'Nx' + str(problem.Nx)
        if problem.Ny is not None:
            s += 'Ny' + str(problem.Ny)

        return s

    def __str__(self):
        return 'Steady NSE'
