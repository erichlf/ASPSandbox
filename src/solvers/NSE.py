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
    def weak_residual(self,U,U_,p,p_,v,q):
        Re = self.Re #Reynolds Number

        alpha = self.alpha #time stepping method
        dt = self.dt

        inviscid = self.options['inviscid']

        #U_(k+alpha)
        U_alpha = (1.0-alpha)*U_ + alpha*U

        #p_(k+alpha)
        p_alpha = (1.0-alpha)*p_ + alpha*p

        #weak form of the equations
        r = (1./dt)*inner(U - U_,v)*dx \
            - p_alpha*div(v)*dx \
            + inner(grad(U_alpha)*U_alpha,v)*dx
        r +=  inviscid/Re*inner(grad(U_alpha),grad(v))*dx
        r += div(U_alpha)*q*dx

        return r

    def stabilization_parameters(self,U_,eta_,h):
        k1  = 1.
        k2  = 0.5
        d1 = k1*(self.dt**(-2) + eta_*eta_*h**(-2))**(-0.5)
        d2 = k2*(self.dt**(-2) + inner(U_,U_)*h**(-2))**(-0.5)

        return d1, d2

    def __str__(self):
          return 'NSE'
