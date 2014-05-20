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
        alpha = self.alpha #time stepping method

        problem = self.problem

        #problem dimensions
        self.x0 = problem.mesh.coordinates()[0][0]
        self.x1 = problem.mesh.coordinates()[-1][0]
        self.y0 = problem.mesh.coordinates()[0][1]
        self.y1 = problem.mesh.coordinates()[-1][1]

        #Scaling Parameters
        g = 1. #Gravity
        lambda0 = self.options['lambda0'] #typical wavelength
        a0 = self.options['a0'] #Typical wave height
        h0 = self.options['h0'] #Typical depth
        sigma = h0/lambda0
        c0 = (h0*g)**(0.5)
        epsilon = a0/h0

        #Physical Parameters
        hd = 1. #Depth [m]
        ad = 0.4 #height of the moving object [m]
        bh = 0.7 #width of the moving object

        #Scaled Parameters
        self.dt = self.dt*c0/lambda0 #time step
        dt = self.dt
        self.T = self.T*c0/lambda0 #Final time
        hd = hd/h0 #depth
        ad = ad/a0 #height of the moving object

        #Definition of the wave_object
        seabed = 'hd'
        zeta0 = seabed
        
        self.zeta = Expression(zeta0, hd=hd, t=self.t0)
        self.zeta_ = Expression(zeta0, hd=hd, t=self.t0)
        self.zeta__ = Expression(zeta0, hd=hd,t=self.t0)
        
        self.z = interpolate(self.zeta, self.Q)
        self.zz_ = interpolate(self.zeta, self.Q)

        zeta_tt = 1./dt**2*(self.zeta - 2*self.zeta_ + self.zeta__)
        zeta_t = 1./dt*(self.zeta - self.zeta_)

        #weak form of the equations

        r = 1./dt*inner(U-U_,v)*dx + epsilon*inner(grad(U)*U,v)*dx \
            - div(v)*eta*dx

        r += sigma**2*1./dt*div(self.zeta*(U-U_))*div(self.zeta*v/2.)*dx \
              - sigma**2*1./dt*div(U-U_)*div(self.zeta**2*v/6.)*dx
        r += sigma**2*zeta_tt/epsilon*div(self.zeta*v/2.)*dx

        r += 1./dt*(eta-eta_)*chi*dx + 1./epsilon*zeta_t*chi*dx
        r -= inner(U,grad(chi))*(epsilon*eta+self.zeta)*dx

        return r

    def wave_object(self, t, dt):
        self.zeta.t = t
        self.zeta_.t = t - dt
        self.zeta__.t = t - 2*dt
        self.zz_.assign(self.z)
        self.z = interpolate(self.zeta, self.Q)

    def __str__(self):
          return 'Peregrine'
