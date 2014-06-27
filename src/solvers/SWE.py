__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"

from solverbase import *

def f(f0,beta): #Coriolis parameter
    return Expression('f0 + beta*x[1]', f0=f0, beta=beta)

class Solver(SolverBase):
    def __init__(self, options):
        SolverBase.__init__(self, options)
        self.Re = options['Re']
        self.H = options['H']
        self.Ro = options['Ro']
        self.Fr = options['Fr']
        self.Th = options['Theta']

        #if we want a linear version then make a coefficient zero for the
        #terms which only occur in the non-linear from of SWE
        if(self.options['linear']):
            self.NonLinear = 0
        else:
            self.NonLinear = 1

        if(self.options['inviscid']):
            self.inviscid = 0
        else:
            self.inviscid = 1


    #strong residual for cG(1)cG(1)
    def strong_residual(self,u,U,eta):
        #get problem parameters
        Re = self.Re #Reynolds number
        H = self.H #Fluid depth
        Ro = self.Ro #Rossby number
        Th = self.Th #average wave height
        Fr = self.Fr #Froude number

        NonLinear = self.NonLinear

        #momentum equation
        R1 = NonLinear*grad(u)*U \
            + 1/Ro*as_vector((-U[1],U[0])) \
            + Fr**(-2)*Th*grad(eta)
        #continuity equation
        R2 = 1/Th*H*div(U)

        return R1, R2

    #weak residual for cG(1)cG(1)
    def weak_residual(self,w,w_,wt,ei_mode):
        (U, eta) = (as_vector((w[0], w[1])), w[2])
        (U_, eta_) = (as_vector((w_[0], w_[1])), w_[2])
        (v, chi) = (as_vector((wt[0], wt[1])), wt[2])

        problem = self.problem

        h = CellSize(self.mesh) #mesh size
        d1, d2 = self.stabilization_parameters(U_,eta_,h) #stabilization parameters

        #set up error indicators
        Z = FunctionSpace(self.mesh, "DG", 0)
        z = TestFunction(Z)

        #get problem parameters
        Re = self.Re #Reynolds number
        H = self.H #Fluid depth
        Ro = self.Ro #Rossby number
        Th = self.Th #average wave height
        Fr = self.Fr #Froude number

        NonLinear = self.NonLinear
        inviscid = self.inviscid

        alpha = self.alpha #time stepping method
        dt = self.dt
        t0 = self.t0

        #U_(k+alpha)
        U_alpha = (1.0-alpha)*U_ + alpha*U

        #p_(k+alpha)
        eta_alpha = (1.0-alpha)*eta_ + alpha*eta

        t = t0 + dt
        #forcing and mass source/sink
        F1_alpha = alpha*problem.F1(t) + (1 - alpha)*problem.F1(t0)
        F2_alpha = alpha*problem.F2(t) + (1 - alpha)*problem.F2(t0)

        #weak form of the equations
        #momentum equation
        r = (1./dt)*inner(U - U_,v)*dx \
            + 1/Ro*(U_alpha[0]*v[1] - U_alpha[1]*v[0])*dx \
            - Fr**(-2)*Th*eta_alpha*div(v)*dx
        r += inviscid/Re*inner(grad(U_alpha),grad(v))*dx
        #add the terms for the non-linear SWE
        r += NonLinear*inner(grad(U_alpha)*U_alpha,v)*dx
        #continuity equation
        r += (1./dt)*(eta - eta_)*chi*dx \
            + H/Th*div(U_alpha)*chi*dx

        r -= inner(F1_alpha,v)*dx + F2_alpha*chi*dx

        #least squares stabilization
        if(not self.options["stabilize"]):
          d1 = 0
          d2 = 0
        if(not ei_mode):
          z = 1.

        R1, R2 = self.strong_residual(U_alpha,U_alpha,eta_alpha)
        Rv1, Rv2 = self.strong_residual(U_alpha,v,chi)
        r += z*(d1*inner(R1 - F1_alpha, Rv1)*dx + d2*(R2 - F2_alpha)*Rv2*dx)

        return r

    def stabilization_parameters(self,U_,eta_,h):
        k1  = (self.Ro*self.Fr**2*self.Th**(-1))/2
        k2  = self.Th/(2*self.H)
        d1 = k1*(self.dt**(-2) + inner(U_,U_)*h**(-1))**(-0.5)
        d2 = k2*(self.dt**(-2) + eta_*eta_*h**(-1))**(-0.5)

        return d1, d2

    def __str__(self):
          return 'SWE'
