__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"

from solverbase import *


class Solver(SolverBase):

    '''
        Solver class for solving the Quasi-Geostrophic equations in
        Steamfunction-Vorticity form.
    '''

    def __init__(self, options):
        SolverBase.__init__(self, options)
        self.Re = options['Re']
        self.Ro = options['Ro']

    def function_space(self, mesh):
        # Define function spaces
        Q = FunctionSpace(mesh, 'CG', self.Pu)
        P = FunctionSpace(mesh, 'CG', self.Pp)
        W = MixedFunctionSpace([Q, P])

        return W

    def weak_residual(self, problem, k, W, w, ww, w_, wt, ei_mode=False):
        (q, psi) = (w[0], w[1])
        (Q, Psi) = (ww[0], ww[1])
        (Q_, Psi_) = (w_[0], w_[1])
        (p, chi) = (wt[0], wt[1])

        h = CellSize(W.mesh())  # mesh size

        # set up error indicators
        Z = FunctionSpace(W.mesh(), "DG", 0)
        z = TestFunction(Z)

        Re = self.Re  # Reynolds Number
        Ro = self.Ro  # Rossby Number

        self.k = Expression('dt', dt=problem.k)
        t0 = problem.t0

        t = t0 + self.k
        # forcing and mass source/sink
        f = problem.F(t)

        # least squares stabilization
        if(not self.options["stabilize"] or ei_mode):
            d1 = 0
            d2 = 0
        if(not ei_mode):
            z = 1.

        # weak form of the equations
        r = z * ((1. / self.k) * (Q - Q_) * p
                 + (1. / Re) * inner(grad(q), grad(p))
                 + self.Jac(psi, q) * p
                 - psi.dx(0) * p) * dx
        r -= z * f * p * dx  # forcing function
        r += z * (q * chi - Ro * inner(grad(psi), grad(chi))) * dx

        return r

    def Jac(self, psi, q):
        return psi.dx(1) * q.dx(0) - psi.dx(0) * q.dx(1)

    def functional(self, mesh, w):

        (q, psi) = (w[0], w[1])

        M = q * dx  # Mean of the vorticity in the whole domain

        return M

    def file_naming(self, n=-1, dual=False):
        if n == -1:
            self._ufile = File(self.s + '_psi.pvd', 'compressed')
            self._pfile = File(self.s + '_q.pvd', 'compressed')
            self._uDualfile = File(self.s + '_psiDual.pvd', 'compressed')
            self._pDualfile = File(self.s + '_qDual.pvd', 'compressed')
            self.meshfile = File(self.s + '_mesh.xml')
        else:
            self._ufile = File(self.s + '_psi%d.pvd' % n, 'compressed')
            self._pfile = File(self.s + '_q%d.pvd' % n, 'compressed')
            self._uDualfile = File(self.s + '_psiDual%d.pvd' % n, 'compressed')
            self._pDualfile = File(self.s + '_qDual%d.pvd' % n, 'compressed')
            self.meshfile = File(self.s + '_mesh%d.xml' % n)

    def __str__(self):
        return 'QGE'
