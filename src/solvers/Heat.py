__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"

from solverbase import *

class InitialConditions(Expression):
    def __init__(self, problem, Q, Psi):
        self.q0, self.psi0 = problem.initial_conditions(Q, Psi)
        self.q0 = project(self.q0,Q)
        self.psi0 = project(self.psi0,Psi)

    def eval(self, value, x):
        value[0] = self.q0(x)
        value[1] = self.psi0(x)

    def value_shape(self):
        return (2,)

class Solver(SolverBase):

    def __init__(self, options):
        SolverBase.__init__(self, options)
        self.a = options['Re']

    def function_space(self, mesh):
        # Define function spaces
        Q = FunctionSpace(mesh, 'CG', self.Pu)

        return Q

    def weak_residual(self, problem, W, w, w_, wt, ei_mode=False):

        h = CellSize(W.mesh()) #mesh size

        #set up error indicators
        Z = FunctionSpace(W.mesh(), "DG", 0)
        z = TestFunction(Z)

        a = self.a #Reynolds Number

        alpha = self.alpha #time stepping method
        k = problem.k
        t0 = problem.t0

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
        r = z*((1./k)*(w - w_)*wt \
            + a*inner(grad(w_alpha),grad(wt)))*dx
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
            self.vizU = plot(w, title='Temperature', rescale=True)
        else :
            self.vizU.plot(w)

    def __str__(self):
          return 'QGE'
