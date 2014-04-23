__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"

from solverbase import *

class Solver(SolverBase):
#    Incremental pressure-correction scheme.

    def __init__(self, options):
        SolverBase.__init__(self, options)

        #parameters
        self.Re = None

    def weak_residual(self,U,U_,eta,eta_,v,chi):
        dt = self.dt #time step
        alpha = self.alpha #time stepping method

        problem = self.problem

        #problem dimensions
        self.x0 = problem.mesh.coordinates()[0][0]
        self.x1 = problem.mesh.coordinates()[-1][0]
        self.y0 = problem.mesh.coordinates()[0][1]
        self.y1 = problem.mesh.coordinates()[-1][1]

        #parameters
        g = 9.8 #Gravity
        lambda0 = 20. #typical wavelength
        a0 = 0.4 #Typical wave height
        h0 = 1. #Typical depth
        sigma = h0/lambda0
        c0 = (h0*g)**(0.5)
        epsilon = a0/h0

        ad = 0.4
        xh = 0.0
        vh = 1.
        bh = 0.7
        hd = 1.
        lamda0=20.

        #self.zeta = self.wave_object(self.t0)
        zeta0 = 'hd - epsilon*ad*exp(-pow(lambda0*(x[0]-xh)-vh*t,2)/(bh*bh))'
        self.zeta = Expression(zeta0, ad=ad, xh=xh, vh=vh, bh=bh, t=self.t0, hd=hd, epsilon=epsilon, lambda0=lambda0, element=self.Q.ufl_element())
        self.zeta_ = Expression(zeta0, ad=ad, xh=xh, vh=vh, bh=bh, t=self.t0, hd=hd, epsilon=epsilon, lambda0=lambda0, element=self.Q.ufl_element())
        self.zeta__ = Expression(zeta0, ad=ad, xh=xh, vh=vh, bh=bh, t=self.t0, hd=hd, epsilon=epsilon, lambda0=lambda0, element=self.Q.ufl_element())

        zeta_tt = 1./dt**2*(self.zeta - 2*self.zeta_ + self.zeta__)
        zeta_t = 1./dt*(self.zeta - self.zeta_)

        #weak form of the equations

        r = 1./dt*inner(U-U_,v)*dx + epsilon*inner(grad(U)*U,v)*dx \
            - div(v)*eta*dx

        r += sigma**2*1./dt*div(self.zeta*(U-U_))*div(self.zeta*v/2.)*dx \
              - sigma**2*1./dt*div(U-U_)*div(self.zeta**2*v/6.)*dx
        r += sigma**2/epsilon*zeta_tt*div(self.zeta*v/2.)*dx

        r += 1./dt*(eta-eta_)*chi*dx + 1./epsilon*zeta_t*chi*dx
        r -= inner(U,grad(chi))*(epsilon*eta+self.zeta)*dx

        return r

    def __str__(self):
          return 'Dimensioinless_Peregrine'
