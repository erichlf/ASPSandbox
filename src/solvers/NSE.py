__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"

from solverbase import *

class Solver(SolverBase):

    def __init__(self, options):
        SolverBase.__init__(self, options)
        self.nu = options['nu']

    #strong residual for cG(1)cG(1)
    def strong_residual(self,u,U,p):
        R1 = div(U)
        R2 = grad(u)*U + grad(p)

        return R1, R2

    #weak residual for cG(1)cG(1)
    def weak_residual(self,U,U_,p,p_,v,q):
        nu = self.nu #Reynolds Number

        theta = self.theta #time stepping method
        dt = self.dt

        #U_(k+theta)
        U_theta = (1.0-theta)*U_ + theta*U

        #p_(k+theta)
        p_theta = (1.0-theta)*p_ + theta*p

        #weak form of the equations
        r = div(U_theta)*q*dx 
        r += (1./dt)*inner(U - U_,v)*dx \
            - p_theta*div(v)*dx \
            + inner(grad(U_theta)*U_theta,v)*dx \
            + nu*inner(grad(U_theta),grad(v))*dx

        return r

    def stabilization_parameters(self,U_,eta_,h):
        k1  = 1.
        k2  = 0.5
        d1 = k1*(self.dt**(-2) + eta_*eta_*h**(-2))**(-0.5)
        d2 = k2*(self.dt**(-2) + inner(U_,U_)*h**(-2))**(-0.5)

        return d1, d2

    def __str__(self):
          return 'NSE'
