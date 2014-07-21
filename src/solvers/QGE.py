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
        self.Re = options['Re']
        self.Ro = options['Ro']

    def function_space(self, mesh):
        # Define function spaces
        Q = FunctionSpace(mesh, 'CG', self.Pu)
        Psi = FunctionSpace(mesh, 'CG', self.Pp)
        W = MixedFunctionSpace([Q, Psi])

        return W

    def weak_residual(self, problem, W, w, w_, wt, ei_mode=False):
        (q, psi) = (w[0], w[1])
        (q_, psi_) = (w_[0], w_[1])
        (p, chi) = (wt[0], wt[1])

        h = CellSize(W.mesh()) #mesh size

        #set up error indicators
        Z = FunctionSpace(W.mesh(), "DG", 0)
        z = TestFunction(Z)

        Re = self.Re #Reynolds Number
        Ro = self.Ro #Rossby Number

        alpha = self.alpha #time stepping method
        k = problem.k
        t0 = problem.t0

        #psi_(k+alpha)
        psi_alpha = (1.0-alpha)*psi_ + alpha*psi
        #q_(k+alpha)
        q_alpha = (1.0-alpha)*q_ + alpha*q

        t = t0 + k
        #forcing and mass source/sink
        f = problem.F(t)

        #least squares stabilization
        if(not self.options["stabilize"] or ei_mode):
          d1 = 0
          d2 = 0
        if(not ei_mode):
          z = 1.

        #weak form of the equations
        r = z*((1./k)*(q - q_)*p \
            + (1./Re)*inner(grad(q_alpha),grad(p)) \
            + self.Jac(psi,q_alpha)*p \
            - psi.dx(0)*p)*dx
        r -= z*f*p*dx #forcing function
        r += z*(q_alpha*chi - Ro*inner(grad(psi),grad(chi)))*dx

        return r

    def Jac(self, psi, q):
        return psi.dx(1)*q.dx(0) - psi.dx(0)*q.dx(1)

    def functional(self, mesh, w):

      (q, psi) = (w[0], w[1])

      M = q*dx # Mean of the x-velocity in the whole domain

      return M

    def __str__(self):
          return 'QGE'
