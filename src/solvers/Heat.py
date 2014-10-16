__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"

from AFES import *
from AFES import Solver as SolverBase


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

        t0 = problem.t0

        kappa = problem.kappa
        rho = problem.rho
        c = problem.rho

        t = t0 + k
        # forcing and mass source/sink
        self.f_ = problem.F(t - k)
        self.f = problem.F(t)

        # weak form of the equations
        r = rho * c * (1. / k) * (U - U_) * v * dx
        if kappa.rank() == 0:
            r += kappa * inner(grad(u), grad(v)) * dx
        else:  # anisotropic case
            r += inner(dot(kappa, grad(u)), grad(v)) * dx

        return r

    def condition(self, ei, m, m_):
        '''
            Adaptive stopping criterion for Galerkin-orthogonal problem (Heat).
            ei - error indicators (non-Galerkin-orthogonal problems)
            m - current functional size (Galerkin-orthogonal problems)
            m_ - previous functional size (Galerkin-orthogonal problems)
        '''
        return abs(m - m_)

    def Save(self, problem, u, dual=False):
        if self.options['save_frequency'] != 0 \
                and (self._timestep - 1) % self.options['save_frequency'] == 0:
            if not dual:
                self._ufile << u
            else:
                self._uDualfile << u

    def file_naming(self, n=-1, dual=False):
        if n == -1:
            self._ufile = File(self.s + '_u.pvd', 'compressed')
            self._uDualfile = File(self.s + '_uDual.pvd', 'compressed')
            self.meshfile = File(self.s + '_mesh.xml')
        else:
            self._ufile = File(self.s + '_u%d.pvd' % n, 'compressed')
            self._uDualfile = File(self.s + '_uDual%d.pvd' % n, 'compressed')
            self.meshfile = File(self.s + '_mesh%d.xml' % n)

    def Plot(self, problem, V, u):

        # Plot velocity and height and wave object
        if self.vizU is None:
            self.vizU = plot(u, title='Temperature', rescale=True, elevate=0.0)
        else:
            self.vizU.plot(u)

    def __str__(self):
        return 'Heat'
