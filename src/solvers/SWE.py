__author__ = "Erich L Foster <efoster@bcamath.org>"
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
        self.Th = options['Th']

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
    def weak_residual(self,U,U_,eta,eta_,v,chi):
        #get problem parameters
        Re = self.Re #Reynolds number
        H = self.H #Fluid depth
        Ro = self.Ro #Rossby number
        Th = self.Th #average wave height
        Fr = self.Fr #Froude number

        NonLinear = self.NonLinear

        alpha = self.alpha #time stepping method
        dt = self.dt

        #U_(k+alpha)
        U_alpha = (1.0-alpha)*U_ + alpha*U

        #p_(k+alpha)
        eta_alpha = (1.0-alpha)*eta_ + alpha*eta

        #weak form of the equations
        #momentum equation
        r = (1./dt)*inner(U - U_,v)*dx \
            + 1/Ro*(U_alpha[0]*v[1] - U_alpha[1]*v[0])*dx \
            - Fr**(-2)*Th*eta_alpha*div(v)*dx \
            + 1/Re*inner(grad(U_alpha),grad(v))*dx
        #add the terms for the non-linear SWE
        r += NonLinear*inner(grad(U_alpha)*U_alpha,v)*dx
        #continuity equation
        r += (1./dt)*(eta - eta_)*chi*dx \
            + H/Th*div(U_alpha)*chi*dx

        return r

    def stabilization_parameters(self,U_,eta_,h):
        k1  = (self.Ro*self.Fr**2*self.Th**(-1))/2
        k2  = self.Th/(2*self.H)
        d1 = k1*(self.dt**(-2) + inner(U_,U_)*h**(-2))**(-0.5)
        d2 = k2*(self.dt**(-2) + eta_*eta_*h**(-2))**(-0.5)

        return d1, d2

    def __str__(self):
          return 'SWE'
