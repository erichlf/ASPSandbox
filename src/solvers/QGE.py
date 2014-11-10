__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"

from AFES import *
from AFES import Solver as SolverBase


class Solver(SolverBase):

    '''
        Solver class for solving the Quasi-Geostrophic equations in
        Steamfunction-Vorticity form.
    '''

    def __init__(self, options):
        try:
            self.Re = options['Re']
        except:
            self.Re = 200
            options['Re'] = self.Re
        try:
            self.Ro = options['Ro']
        except:
            self.Ro = 0.0016
            options['Ro'] = self.Ro

        SolverBase.__init__(self, options)

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

        Re = self.Re  # Reynolds Number
        Ro = self.Ro  # Rossby Number

        t0 = problem.t0

        t = t0 + k
        # forcing and mass source/sink
        f = problem.F(t)

        # weak form of the equations
        r = ((1. / k) * (Q - Q_) * p
             + (1. / Re) * inner(grad(q), grad(p))
             + self.Jac(psi, q) * p
             - psi.dx(0) * p) * dx
        r -= f * p * dx  # forcing function
        r += (q * chi - Ro * inner(grad(psi), grad(chi))) * dx

        return r

    def Jac(self, psi, q):
        return psi.dx(1) * q.dx(0) - psi.dx(0) * q.dx(1)

    def file_naming(self, n=-1, dual=False):
        if n == -1:
            self._ufile = File(self.s + '_q.pvd', 'compressed')
            self._pfile = File(self.s + '_psi.pvd', 'compressed')
            self._uDualfile = File(self.s + '_qDual.pvd', 'compressed')
            self._pDualfile = File(self.s + '_psiDual.pvd', 'compressed')
            self.meshfile = File(self.s + '_mesh.xml')
        else:
            self._ufile = File(self.s + '_q%d.pvd' % n, 'compressed')
            self._pfile = File(self.s + '_psi%d.pvd' % n, 'compressed')
            self._uDualfile = File(self.s + '_qDual%d.pvd' % n, 'compressed')
            self._pDualfile = File(self.s + '_psiDual%d.pvd' % n, 'compressed')
            self.meshfile = File(self.s + '_mesh%d.xml' % n)

    def __str__(self):
        return 'QGE'
