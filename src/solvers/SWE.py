__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"

from solverbase import *

def f(f0,beta): #Coriolis parameter
    return Expression('f0 + beta*x[1]', f0=f0, beta=beta)

class Solver(SolverBase):
    def __init__(self, options):
        SolverBase.__init__(self, options)

    #strong residual for cG(1)cG(1)
    def strong_residual(self,u,U,eta):
        H = 1#self.H
        f0 = self.f0
        beta = self.beta
        g = 1.#self.g
        nu = self.nu
        rho = self.rho

        NonLinear = self.NonLinear

        R1 = H*div(U)
        R2 = NonLinear*grad(u)*U \
            + f(f0,beta)*as_vector((-U[1],U[0])) \
            + g*grad(eta)

        return R1, R2

    #weak residual for cG(1)cG(1)
    def weak_residual(self,U,U_,eta,eta_,v,chi):
        H = self.H
        f0 = self.f0
        beta = self.beta
        g = self.g
        nu = self.nu
        rho = self.rho

        NonLinear = self.NonLinear

        theta = self.options['theta'] #time stepping method
        dt = self.options["dt"]

        #U_(k+theta)
        U_theta = (1.0-theta)*U_ + theta*U

        #p_(k+theta)
        eta_theta = (1.0-theta)*eta_ + theta*eta

        #weak form of the equations
        r = (1./dt)*(eta - eta_)*chi*dx \
            + H*div(U_theta)*chi*dx 
        r += (1./dt)*inner(U - U_,v)*dx \
            + f(f0,beta)*(U_theta[0]*v[1] - U_theta[1]*v[0])*dx \
            - g*eta_theta*div(v)*dx
        #add the terms for the non-linear SWE
        r += NonLinear*(inner(grad(U_theta)*U_theta,v)*dx \
            + nu*inner(grad(U_theta),grad(v))*dx)

        return r

    def __str__(self):
          return 'SWE'
