__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"

from solverbase import *

def f(f0,beta): #Coriolis parameter
    return Expression('f0 + beta*x[1]', f0=f0, beta=beta)

class Solver(SolverBase):
    def __init__(self, options):
        SolverBase.__init__(self, options)
        self.nu = options['nu']
        self.H = options['H']
        self.f0 = options['f0']
        self.beta = options['beta']
        self.g = options['g']

    #strong residual for cG(1)cG(1)
    def strong_residual(self,u,U,eta):
        #get problem parameters
        nu = self.nu
        H = self.H
        f0 = self.f0
        beta = self.beta
        g = self.g

        NonLinear = self.NonLinear

        #assume H is constant, so it doesn't show up in
        #non-dimensional form
        R1 = NonLinear*grad(u)*U \
            + f(f0,beta)*as_vector((-U[1],U[0])) \
            + g*grad(eta)
        R2 = H*div(U)

        return R1, R2

    #weak residual for cG(1)cG(1)
    def weak_residual(self,U,U_,eta,eta_,v,chi):
        #get problem parameters
        nu = self.nu
        H = self.H
        f0 = self.f0
        beta = self.beta
        g = self.g

        NonLinear = self.NonLinear

        theta = self.theta #time stepping method
        dt = self.dt

        #U_(k+theta)
        U_theta = (1.0-theta)*U_ + theta*U

        #p_(k+theta)
        eta_theta = (1.0-theta)*eta_ + theta*eta

        #weak form of the equations
        r = (1./dt)*inner(U - U_,v)*dx \
            + f(f0,beta)*(U_theta[0]*v[1] - U_theta[1]*v[0])*dx \
            - g*eta_theta*div(v)*dx \
            + nu*inner(grad(U_theta),grad(v))*dx
        #add the terms for the non-linear SWE
        r += NonLinear*inner(grad(U_theta)*U_theta,v)*dx
        r += (1./dt)*(eta - eta_)*chi*dx \
            + H*div(U_theta)*chi*dx

        return r

    def stabilization_parameters(self,U_,eta_,h):
        k1  = 1./(2*self.g)
        k2  = 1./(2*self.H)
        d1 = k1*(self.dt**(-2) + eta_*eta_*h**(-2))**(-0.5)
        d2 = k2*(self.dt**(-2) + inner(U_,U_)*h**(-2))**(-0.5)

        return d1, d2

    def __str__(self):
          return 'SWE'
