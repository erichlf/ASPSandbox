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
        SolverBase.__init__(self, options)

    def function_space(self, mesh):
        # Define function spaces
        V = FunctionSpace(mesh, 'CG', self.Pu)

        return V

    # strong residual for cG(1)cG(1)
    def strong_residual(self, alpha, beta, u):

        R = alpha * u + dot(beta, grad(u))

        return R

    # weak residual for cG(1)cG(1)
    def weak_residual(self, problem, k, V, u, U, U_, v, ei_mode=False):
        h = CellSize(V.mesh())  # mesh size

        alpha = problem.alpha  # reaction coefficient
        kappa = problem.kappa
        beta = problem.beta  # velocity

        d = self.stabilization_parameters(U_, h, beta, kappa)

        t0 = problem.t0

        t = t0 + k
        # forcing and mass source/sink
        F = problem.F(t)

        # least squares stabilization
        if ei_mode:
            d = Constant(0)

        # weak form of the equations
        r = ((1. / k) * inner(U - U_, v)
             + alpha * inner(u, v)) * dx

        if kappa.rank() == 0:
            r += kappa * inner(grad(u), grad(v)) * dx
        else:  # anisotropic case
            r += inner(dot(kappa, grad(u)), grad(v)) * dx

        r += inner(dot(beta, grad(u)), v) * dx

        # forcing function
        r -= inner(F, v) * dx

        r += d * inner(self.strong_residual(alpha, beta, u) - F,
                       self.strong_residual(alpha, beta, v)) * dx

        return r

    def stabilization_parameters(self, u, h, beta, kappa):
        d = conditional(ge(h, kappa), h/sqrt(dot(beta, beta)),
                h**2/sqrt(dot(beta, beta)))

        return d

    def Save(self, problem, u, dual=False):
        if self.saveFrequency != 0 \
                and (self._timestep - 1) % self.saveFrequency == 0:
            if not dual:
                u.rename("primal", "ADR")
                self._ufile << u
            else:
                u.rename("dual", "ADR")
                self._uDualfile << u

    def file_naming(self, problem, n=-1, opt=False):
        s = 'results/' + self.prefix(problem) + self.suffix(problem)

        if n == -1:
            if opt:
                self._ufile = File(s + '_uOpt.pvd', 'compressed')
            else:
                self._ufile = File(s + '_u.pvd', 'compressed')
            self._uDualfile = File(s + '_uDual.pvd', 'compressed')
            self.meshfile = File(s + '_mesh.xml')
        else:
            self.eifile = File(s + '_ei.pvd', 'compressed')
            self._ufile = File(s + '_u%d.pvd' % n, 'compressed')
            self._uDualfile = File(s + '_uDual%d.pvd' % n, 'compressed')
            self.meshfile = File(s + '_mesh%d.xml' % n)

    def suffix(self, problem):
        import numpy as np
        '''
            Obtains the run specific data for file naming, e.g. Nx, k, etc.
        '''
        s = ''

        alpha = float(problem.alpha)
        beta = np.zeros(2)
        problem.beta.eval(beta, np.zeros(3))
        beta = np.dot(beta, beta)**0.5
        if problem.kappa.rank() == 0:
            kappa = float(problem.kappa)
        else:
            kappa = np.zeros(2)
            problem.kappa.eval(kappa, np.zeros(3))
            kappa = np.dot(kappa, kappa)**0.5

        # Return file suffix for output files
        s = 'alpha%.3Gbeta%.3Gkappa%.3G' % (alpha, beta, kappa)

        s += 'T%G' % (problem.T)
        if problem.Nx is not None:
            s += 'Nx' + str(problem.Nx)
        if problem.Ny is not None:
            s += 'Ny' + str(problem.Ny)
        if problem.mesh.topology().dim() > 2 and problem.Nz is not None:
            s += 'Nz' + str(problem.Nz)

        s += 'K' + str(problem.k)

        return s

    def Plot(self, problem, V, u):

        # Plot velocity and height and wave object
        if self.vizU is None:
            self.vizU = plot(u, title='Concentration', rescale=True,
                             elevate=0.0)
        else:
            self.vizU.plot(u)

    def __str__(self):
        return 'ADR'
