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
        #get problem parameters
        Fr = self.Fr
        Th = self.Th
        Ro = self.Ro
        Re = self.Re

        NonLinear = self.NonLinear

        #assume H is constant, so it doesn't show up in
        #non-dimensional form
        R1 = 1./Th*div(U)
        R2 = NonLinear*grad(u)*U \
            + 1./Ro*as_vector((-U[1],U[0])) \
            + Fr**(-2.)*Th*grad(eta)

        return R1, R2

    #weak residual for cG(1)cG(1)
    def weak_residual(self,U,U_,eta,eta_,v,chi):
        #get problem parameters
        Fr = self.Fr
        Th = self.Th
        Ro = self.Ro
        Re = self.Re

        NonLinear = self.NonLinear

        theta = self.options['theta'] #time stepping method
        dt = self.options["dt"]

        #U_(k+theta)
        U_theta = (1.0-theta)*U_ + theta*U

        #p_(k+theta)
        eta_theta = (1.0-theta)*eta_ + theta*eta

        #weak form of the equations
        #assume H is constant, so it doesn't show up in
        #non-dimensional form
        r = (1./dt)*(eta - eta_)*chi*dx \
            + 1./Th*div(U_theta)*chi*dx 
        r += (1./dt)*inner(U - U_,v)*dx \
            + 1./Ro*(U_theta[0]*v[1] - U_theta[1]*v[0])*dx \
            - Fr**(-2.)*Th*eta_theta*div(v)*dx
        #add the terms for the non-linear SWE
        r += NonLinear*(inner(grad(U_theta)*U_theta,v)*dx \
            + 1./Re*inner(grad(U_theta),grad(v))*dx)

        return r

    def __str__(self):
          return 'SWE'
