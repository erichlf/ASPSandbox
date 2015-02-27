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
        self.steady_state = True

    def function_space(self, mesh):
        # Define function spaces
        V = FunctionSpace(mesh, 'CG', 1)

        return V

    def weak_residual(self, problem, W, u, v, ei_mode=False):
        # problem parameters
        kappa = problem.kappa

        self.f = problem.F()  # forcing and mass source/sink

        # weak form of the equations
        if kappa.rank() == 0:
            r = kappa * inner(grad(u), grad(v)) * dx
        else:  # anisotropic case
            r = inner(dot(kappa, grad(u)), grad(v)) * dx

        r -= self.f * v * dx

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
        if self.saveFrequency != 0 \
                and (self._timestep - 1) % self.saveFrequency == 0:
            if not dual:
                self._ufile << u
            else:
                self._uDualfile << u

    def suffix(self, problem):
        import numpy as np
        '''
            Obtains the run specific data for file naming, e.g. Nx, k, etc.
        '''
        s = ''

        if problem.kappa.rank() == 0:
            kappa = float(problem.kappa)
        else:
            kappa = np.zeros(2)
            problem.kappa.eval(kappa, np.zeros(3))
            kappa = np.dot(kappa, kappa)**0.5

        # Return file suffix for output files
        s = 'kappa%.3G' % (kappa)

        if problem.Nx is not None:
            s += 'Nx' + str(problem.Nx)
        if problem.Ny is not None:
            s += 'Ny' + str(problem.Ny)
        if problem.mesh.topology().dim() > 2 and problem.Nz is not None:
            s += 'Nz' + str(problem.Nz)

        return s

    def file_naming(self, problem, n=-1, opt=False):
        s = self.dir + self.prefix(problem) + self.suffix(problem)

        if n == -1:
            if opt:
                self._ufile = File(s + '_uOpt.pvd', 'compressed')
            else:
                self._ufile = File(s + '_u.pvd', 'compressed')
            self.eifile = File(s + '_ei.pvd', 'compressed')
            self._uDualfile = File(s + '_uDual.pvd', 'compressed')
            self.meshfile = File(s + '_mesh.xml')
            self.optfile = File(s + '_control.pvd', 'compressed')
        else:
            if self.eifile is None:  # error indicators
                self.eifile = File(s + '_ei.pvd', 'compressed')
            self._ufile = File(s + '_u%d.pvd' % n, 'compressed')
            self._uDualfile = File(s + '_uDual%d.pvd' % n, 'compressed')
            self.meshfile = File(s + '_mesh%d.xml' % n)

    def Plot(self, problem, V, u):

        # Plot velocity and height and wave object
        if self.vizU is None:
            self.vizU = plot(u, title='Temperature', rescale=False)
        else:
            self.vizU.plot(u)

    def __str__(self):
        return 'Poisson'
