__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"

from solverbase import *

class Solver(SolverBase):
#    Incremental pressure-correction scheme.

    def __init__(self, options):
        SolverBase.__init__(self, options)

    #strong residual for cG(1)cG(1)
    def strong_residual(self,u,U,p):
        nu = self.nu
        rho = self.rho #density

        R1 = div(U)
        R2 = grad(u)*U + grad(p)

        return R1, R2

    #weak residual for cG(1)cG(1)
    def weak_residual(self,U,U_,p,p_,v,q):
        nu = self.nu
        rho = self.rho #density

        theta = self.options['theta'] #time stepping method
        dt = self.options["dt"]

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

    def __str__(self):
          return 'NSE'
