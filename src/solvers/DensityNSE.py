__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"

from solverbase import *

class Solver(SolverBase):

    def __init__(self, options):
        SolverBase.__init__(self, options)
        self.Re = options['Re']

    # Define function spaces
    def function_space(self, mesh):
        V = VectorFunctionSpace(mesh, 'CG', self.Pu)
        R = FunctionSpace(mesh, 'CG', self.Pp)
        Q = FunctionSpace(mesh, 'CG', self.Pp)

        W = MixedFunctionSpace([V, R, Q])

        return W

    #strong residual for cG(1)cG(1)
    def strong_residual(self,u,U,rho,Rho,p):
        R1 = rho*grad(U)*u + grad(p)
        R2 = inner(U,grad(rho))
        R3 = div(U)

        return R1, R2, R3

    #weak residual for cG(1)cG(1)
    def weak_residual(self,problem,W,w,w_,wt,ei_mode=False):
        #rho = 1/rho' - 1
        (U, rho, p) = (as_vector((w[0], w[1])), w[2], w[3])
        (U_, rho_, p_) = (as_vector((w_[0], w_[1])), w_[2], w_[3])
        (v, r, q) = (as_vector((wt[0], wt[1])), wt[2], wt[3])

        h = CellSize(W.mesh()) #mesh size
        d1, d2, d3 = self.stabilization_parameters(U_, rho_, p_, h) #stabilization parameters
        d = 0.1*h**(3./2.) #stabilization parameter

        #set up error indicators
        Z = FunctionSpace(W.mesh(), "DG", 0)
        z = TestFunction(Z)

        Re = self.Re #Reynolds Number

        alpha = self.alpha #time stepping method
        k = problem.k
        t0 = problem.t0

        #U_(k+alpha)
        U_alpha = (1.0-alpha)*U_ + alpha*U
        rho_alpha = (1.0-alpha)*rho_ + alpha*rho

        t = t0 + k
        #forcing and mass source/sink
        f = problem.F1(t)

        #least squares stabilization
        if(not self.options["stabilize"] or ei_mode):
          d1 = 0
          d2 = 0
          d3 = 0
          d = 0
        if(not ei_mode):
          z = 1.

        #weak form of the equations
        r = z*((1./k)*(rho - rho_) + inner(U_alpha,grad(rho_alpha)))*r*dx #mass equation
        r += z*(rho_alpha*((1./k)*inner(U - U_,v) \
            + inner(grad(U_alpha)*U_alpha,v)) \
            - p*div(v) \
            + 1./Re*inner(grad(U_alpha),grad(v)) \
            + rho_alpha*dot(f,v))*dx #momentum equation
        r += z*div(U_alpha)*q*dx #continuity equation

        R1, R2, R3 = self.strong_residual(U_alpha,U_alpha,rho_alpha,rho_alpha,p)
        Rv1, Rv2, Rv3 = self.strong_residual(U_alpha,v,rho_alpha,r,q)
        r += z*(d1*inner(R1, Rv1) + d2*R2*Rv2 + d3*R3*Rv3)*dx
        r += z*d*(inner(grad(U_alpha),grad(v)))*dx
        #r += z*d*(inner(grad(rho_alpha),grad(r)))*dx

        return r

    def functional(self,mesh,w):

      (u, p) = (as_vector((w[0], w[1])), w[2])

      M = u[0]*dx # Mean of the x-velocity in the whole domain

      return M

    def stabilization_parameters(self,U,rho,p,h):
        K1  = 1.
        K2  = 1.
        K3  = 0.5
        d1 = K1*(self.k**(-2) + inner(U,U)*h**(-2))**(-0.5)
        d2 = K2*(self.k**(-2) + rho*rho*h**(-2))**(-0.5)
        d3 = K3*h**(-1)

        return d1, d2, d3

    def project(self,f1,f2,f3,W):
        #This function will project an expressions into W
        e1, e2, e3 = TestFunctions(W)

        w = TestFunction(W)
        u = Function(W)

        A = inner(u,w)*dx - inner(f1,e1)*dx - f2*e2*dx - f3*e3*dx

        solve(A == 0, u, annotate=False)

        return u

    def __str__(self):
          return 'DensityNSE'
