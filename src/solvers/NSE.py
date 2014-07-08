__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"

from solverbase import *

class Solver(SolverBase):

    def __init__(self, options):
        SolverBase.__init__(self, options)
        self.Re = options['Re']

    #strong residual for cG(1)cG(1)
    def strong_residual(self,u,U,p):
        R1 = grad(u)*U + grad(p)
        R2 = div(U)

        return R1, R2

    #weak residual for cG(1)cG(1)
    def weak_residual(self,W,w,w_,wt,ei_mode=False):
        (U, p) = (as_vector((w[0], w[1])), w[2])
        (U_, p_) = (as_vector((w_[0], w_[1])), w_[2])
        (v, q) = (as_vector((wt[0], wt[1])), wt[2])

        h = CellSize(W.mesh()) #mesh size
        d1, d2 = self.stabilization_parameters(U_,p_,h) #stabilization parameters

        #set up error indicators
        Z = FunctionSpace(W.mesh(), "DG", 0)
        z = TestFunction(Z)

        Re = self.Re #Reynolds Number

        problem = self.problem

        alpha = self.alpha #time stepping method
        k = problem.k
        t0 = problem.t0

        inviscid = self.options['inviscid']

        #U_(k+alpha)
        U_alpha = (1.0-alpha)*U_ + alpha*U

        #p_(k+alpha)
        p_alpha = (1.0-alpha)*p_ + alpha*p

        t = t0 + k
        #forcing and mass source/sink
        F1_alpha = alpha*problem.F1(t) + (1 - alpha)*problem.F1(t0)
        F2_alpha = alpha*problem.F2(t) + (1 - alpha)*problem.F2(t0)

        #least squares stabilization
        if(not self.options["stabilize"] or ei_mode):
          d1 = 0
          d2 = 0
        if(not ei_mode):
          z = 1.

        #weak form of the equations
        r = z*((1./k)*inner(U - U_,v) \
            - p_alpha*div(v) \
            + inner(grad(U_alpha)*U_alpha,v))*dx
        r += z*inviscid/Re*inner(grad(U_alpha),grad(v))*dx
        r += z*div(U_alpha)*q*dx

        r -= z*(inner(F1_alpha,v) + F2_alpha*q)*dx

        R1, R2 = self.strong_residual(U_alpha,U_alpha,p_alpha)
        Rv1, Rv2 = self.strong_residual(U_alpha,v,q)
        r += z*(d1*inner(R1 - F1_alpha, Rv1) + d2*(R2 - F2_alpha)*Rv2)*dx

        return r

    def functional(self,mesh,u,eta):

      M = u[0]*dx # Mean of the x-velocity in the whole domain

      return M

    def stabilization_parameters(self,U_,eta_,h):
        K1  = 1.
        K2  = 0.5
        d1 = K1*(self.k**(-2) + eta_*eta_*h**(-2))**(-0.5)
        d2 = K2*(self.k**(-2) + inner(U_,U_)*h**(-2))**(-0.5)

        return d1, d2

    def __str__(self):
          return 'NSE'
