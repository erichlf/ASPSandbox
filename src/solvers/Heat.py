__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"

from solverbase import *

class Solver(SolverBase):

    def __init__(self, options):
        SolverBase.__init__(self, options)

    def function_space(self, mesh):
        # Define function spaces
        Q = FunctionSpace(mesh, 'CG', self.Pu)

        return Q

    def weak_residual(self, problem, W, w, w_, wt, ei_mode=False):

        h = CellSize(W.mesh()) #mesh size

        #set up error indicators
        Z = FunctionSpace(W.mesh(), "DG", 0)
        z = TestFunction(Z)

        alpha = self.alpha #time stepping method
        k = problem.k
        t0 = problem.t0

        kappa = problem.kappa
        rho = problem.rho
        c = problem.rho

        #u_(k+alpha)
        w_alpha = (1.0-alpha)*w_ + alpha*w

        t = t0 + k
        #forcing and mass source/sink
        self.f_ = problem.F1(t)
        self.f = problem.F1(t)

        f_alpha = (1.0-alpha)*self.f_ + alpha*self.f

        if(not ei_mode):
          z = 1.

        #weak form of the equations
        r = z*(rho*c*(1./k)*(w - w_)*wt \
            + kappa*inner(grad(w_alpha),grad(wt)))*dx
        #r -= z*f_alpha*wt*dx #forcing function

        return r

    def functional(self, mesh, w):

      M = w*dx # Mean of the vorticity in the whole domain

      return M

    def Save(self, problem, w, dual=False):
        k = self.k
        Nx = self.options['Nx']
        Ny = self.options['Ny']
        if (self._timestep - 1) % self.options['save_frequency'] == 0:
            # Create files for saving
            self.file_naming(problem, k, Nx, Ny, dual=False)

            if not dual:
                self._ufile << w
            else:
                self._uDualfile << w

    def file_naming(self, n=-1, dual=False):
        if n==-1:
            self._ufile = File(self.s + '_u.pvd', 'compressed')
            self._uDualfile = File(self.s + '_uDual.pvd', 'compressed')
        else:
            self._ufile = File(self.s + '_u%d.pvd' % n, 'compressed')
            self._uDualfile = File(self.s + '_uDual%d.pvd' % n, 'compressed')

    def Plot(self, problem, W, w):

        # Plot velocity and height and wave object
        if self.vizU is None:
            self.vizU = plot(w, title='Temperature', rescale=True, elevate=0.0)
        else :
            self.vizU.plot(w)

    def __str__(self):
          return 'QGE'
