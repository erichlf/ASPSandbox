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

        d = self.stabilization_parameters(U_, k, h, kappa)

        t0 = problem.t0

        t = t0 + k
        # forcing and mass source/sink
        F = problem.F(t)

        # least squares stabilization
        if ei_mode:
            d = 0

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

    def stabilization_parameters(self, u, k, h, kappa):
        K = 0.5 / sqrt(1.0**2 + 1.61**2)
        d = conditional(le(h, kappa), K*h**2, K*h)

        d = K * h

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
            self.vizU = plot(u, title='Concentration', rescale=True,
                             elevate=0.0)
        else:
            self.vizU.plot(u)

    def __str__(self):
        return 'ADR'
