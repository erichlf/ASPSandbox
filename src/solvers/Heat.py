__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"

from solverbase import *
import numpy as np

class Solver(SolverBase):
    '''
        Solver class for solving the Heat equation.
    '''

    def __init__(self, options):
        SolverBase.__init__(self, options)

    def function_space(self, mesh):
        # Define function spaces
        Q = FunctionSpace(mesh, 'CG', self.Pu)

        return Q

    def weak_residual(self, problem, k, W, w, w_, wt, ei_mode=False):

        h = CellSize(W.mesh()) #mesh size

        #set up error indicators
        Z = FunctionSpace(W.mesh(), "DG", 0)
        z = TestFunction(Z)

        alpha = self.alpha #time stepping method
        t0 = problem.t0

        kappa = problem.kappa
        rho = problem.rho
        c = problem.rho

        #u_(k+alpha)
        w_alpha = (1.0-alpha)*w_ + alpha*w

        t = t0 + self.k
        #forcing and mass source/sink
        self.f_ = problem.F1(t)
        self.f = problem.F1(t)

        if(not ei_mode):
          z = 1.

        #weak form of the equations
        r = z*rho*c*(1./k)*(w - w_)*wt*dx
        if kappa.rank() == 0:
            r += z*kappa*inner(grad(w_alpha),grad(wt))*dx
        else: #anisotropic case
            r += z*inner(dot(kappa,grad(w_alpha)),grad(wt))*dx

        return r

    def functional(self, mesh, w):

      M = w*dx # Mean of the vorticity in the whole domain

      return M

    def condition(self, ei, m, m_):
        '''
            Adaptive stopping criterion for Galerkin-orthogonal problem (Heat).
            ei - error indicators (non-Galerkin-orthogonal problems)
            m - current functional size (Galerkin-orthogonal problems)
            m_ - previous functional size (Galerkin-orthogonal problems)
        '''
        return abs(m - m_)

    def Save(self, problem, w, dual=False):
        if self.options['save_frequency'] !=0 and (self._timestep - 1) % self.options['save_frequency'] == 0:
            if not dual:
                self._ufile << w
            else:
                self._uDualfile << w

    def file_naming(self, n=-1, dual=False):
        if n==-1:
            self._ufile = File(self.s + '_u.pvd', 'compressed')
            self._uDualfile = File(self.s + '_uDual.pvd', 'compressed')
            self.meshfile = File(self.s + '_mesh.xml')
        else:
            self._ufile = File(self.s + '_u%d.pvd' % n, 'compressed')
            self._uDualfile = File(self.s + '_uDual%d.pvd' % n, 'compressed')
            self.meshfile = File(self.s + '_mesh%d.xml' % n)

    def Plot(self, problem, W, w):

        # Plot velocity and height and wave object
        if self.vizU is None:
            self.vizU = plot(w, title='Temperature', rescale=True, elevate=0.0)
        else :
            self.vizU.plot(w)

    def __str__(self):
          return 'Heat'
