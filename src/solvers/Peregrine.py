__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"

from solverbase import *

class Solver(SolverBase):
#    Incremental pressure-correction scheme.

    def __init__(self, options):
        SolverBase.__init__(self, options)
        
        #Parameters
        self.Re = None

    def strong_residual(self,u,U,eta):  
        #Parameters
        g = 9.8 #Gravity
        lambda0 = options['lambda0'] #typical wavelength
        a0 = options['a0'] #Typical wave height
        h0 = options['h0'] #Typical depth
        sigma = h0/lambda0
        c0 = (h0*g)**(0.5)
        epsilon = a0/h0
        
        zeta_tt = 1./self.dt**2*(self.zeta - 2*self.zeta_ + self.zeta__)
        zeta_t = 1./self.dt*(self.zeta - self.zeta_)

        #strong form for stabilization
        R1 = epsilon*grad(u)*U + grad(eta)
        z1 = -sigma**2/2.*epsilon*self.zeta*grad(zeta_tt)

        R2 = div(u)*(epsilon*eta) + inner(u,grad(epsilon*eta)) \
            + div(U)*epsilon*self.zeta + inner(U,grad(epsilon*self.zeta))
        z2 = zeta_t

        return R1, R2, z1, z2

    def weak_residual(self,U,U_,eta,eta_,v,chi):
        alpha = self.alpha #time stepping method
        
        #Parameters
        g = 9.8 #Gravity
        lambda0 = self.options['lambda0'] #typical wavelength
        a0 = self.options['a0'] #Typical wave height
        h0 = self.options['h0'] #Typical depth
        sigma = h0/lambda0
        c0 = (h0*g)**(0.5)
        epsilon = a0/h0
        
        problem = self.problem

        #problem dimensions
        self.x0 = problem.mesh.coordinates()[0][0]
        self.x1 = problem.mesh.coordinates()[-1][0]
        self.y0 = problem.mesh.coordinates()[0][1]
        self.y1 = problem.mesh.coordinates()[-1][1]

        #Scaled Parameters
        self.dt = self.dt*c0/lambda0 #time step
        dt = self.dt
        self.T = self.T*c0/lambda0 #Final time
        
        D = Expression(problem.D, hd=problem.hd, hb=problem.hb,  lambda0=lambda0, element=self.Q.ufl_element()) 
        self.zeta = Expression(problem.zeta0, ad=problem.ad, c0=c0, t0=self.t0, t=self.t0,\
                    hd=problem.hd, lambda0=lambda0, vmax=problem.vmax, element=self.Q.ufl_element())
        self.zeta_ = Expression(problem.zeta0, ad=problem.ad, c0=c0, t0=self.t0, t=self.t0,\
                    hd=problem.hd, lambda0=lambda0, vmax=problem.vmax, element=self.Q.ufl_element())
        self.zeta__ = Expression(problem.zeta0, ad=problem.ad, c0=c0, t0=self.t0, t=self.t0,\
                    hd=problem.hd, lambda0=lambda0, vmax=problem.vmax, element=self.Q.ufl_element())   
        self.bottom = Expression(problem.D + ' + epsilon*(' + problem.zeta0 +')', epsilon=epsilon,\
                    hb=problem.hb, ad=problem.ad, c0=c0, t0=self.t0, t=self.t0, hd=problem.hd,\
                    lambda0=lambda0, vmax=problem.vmax, element=self.Q.ufl_element())
        self.h = interpolate(self.bottom, self.Q)
        self.h_ = interpolate(self.bottom, self.Q)
       
        zeta_tt = 1./dt**2*(self.zeta - 2*self.zeta_ + self.zeta__)
        zeta_t = 1./dt*(self.zeta - self.zeta_)
        
        #Define function to stabilize the wake of the object
        self.time_stabilize = Expression(problem.filtre, vmax=problem.vmax, t=self.t0, lambda0=lambda0, c0=c0)
        
        #Time stepping method
        U_alpha = (1. - alpha)*U_ + alpha*U
        eta_alpha = (1. - alpha)*eta_ + alpha*eta
        zeta_alpha = (1. - alpha)*self.zeta_ + alpha*self.zeta

        #weak form of the equations
        r = 1./dt*inner(U-U_,v)*dx + epsilon*inner(grad(U_alpha)*U_alpha,v)*dx \
            - div(v)*eta_alpha*dx

        r += sigma**2*1./dt*div((D + epsilon*zeta_alpha)*(U-U_))*div((D + epsilon*zeta_alpha)*v/2.)*dx \
              - sigma**2*1./dt*div(U-U_)*div((D + epsilon*zeta_alpha)**2*v/6.)*dx
        r += sigma**2*zeta_tt*div((D + epsilon*zeta_alpha)*v/2.)*dx

        r += 1./dt*(eta-eta_)*chi*dx + zeta_t*chi*dx
        r -= inner(U_alpha,grad(chi))*(epsilon*eta_alpha + D + epsilon*zeta_alpha)*dx

        return r

    def stabilization_parameters(self,U_,eta_,h):
        k1  = 2.
        k2  = 2.
        d1 = k1*(self.dt**(-2) + inner(U_,U_)*h**(-2))**(-0.5)
        d2 = k2*(self.dt**(-2) + eta_*eta_*h**(-2))**(-0.5)

        return d1, d2
    
    def wave_object(self, t, dt): 
        self.zeta.t = t
        self.zeta_.t = t - dt
        self.zeta__.t = t - 2*dt
        self.h_.assign(self.h)
        self.bottom.t = t
        self.h = interpolate(self.bottom, self.Q)
        self.time_stabilize.t = t

    def __str__(self):
          return 'Peregrine'
