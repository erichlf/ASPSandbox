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

    def strong_residual(self,u,U,eta):
        #Scaling Parameters
        g = 9.8 #Gravity
        lambda0 = self.options['lambda0'] #typical wavelength
        a0 = self.options['a0'] #Typical wave height
        h0 = self.options['h0'] #Typical depth
        sigma = h0/lambda0
        c0 = (h0*g)**(0.5)
        epsilon = a0/h0

        #strong form for stabilization
        R1 = epsilon*grad(u)*U + grad(eta)
        z1 = 0.

        R2 = div(u)*(epsilon*eta) + inner(u,grad(epsilon*eta)) \
            + div(U)*epsilon*self.zeta + inner(U,grad(epsilon*self.zeta))
        z2 = 0.

        return R1, R2, z1, z2

    def weak_residual(self,U,U_,eta,eta_,v,chi):
        alpha = self.alpha #time stepping method

        problem = self.problem

        #problem dimensions
        self.x0 = problem.mesh.coordinates()[0][0]
        self.x1 = problem.mesh.coordinates()[-1][0]
        self.y0 = problem.mesh.coordinates()[0][1]
        self.y1 = problem.mesh.coordinates()[-1][1]

        #Scaling Parameters
        g = 9.8 #Gravity
        lambda0 = self.options['lambda0'] #typical wavelength
        a0 = self.options['a0'] #Typical wave height
        h0 = self.options['h0'] #Typical depth
        sigma = h0/lambda0
        c0 = (h0*g)**(0.5)
        epsilon = a0/h0

        #Physical Parameters
        hd = 1. #Depth [m]
        ad = 0.4 #height of the moving object [m]

        #Scaled Parameters
        self.dt = self.dt*c0/lambda0 #time step
        dt = self.dt
        self.T = self.T*c0/lambda0 #Final time
        hd = hd/h0 #depth
        ad = ad/a0 #height of the moving object

        #Definition of the wave_object
        vmax = ((hd*h0+ad*a0)*g)**(0.5) #Max Speed of the moving object [m.s^(-1)]
        seabed = 'hd - 0.5/10.*(x[1]>4./lambda0 ? 1. : 0.)*(lambda0*x[1]-4.) + 0.5/10.*(x[1]<(-4./lambda0) ? 1. : 0.)*(lambda0*x[1]+4.)'
        movingObject = ' - (x[1]<3/lambda0 ? 1. : 0.)*(x[1]>0 ? 1. : 0.)*(lambda0*x[0]>-6 ? 1. : 0.)*ad*0.5*0.5*(1. - tanh(0.5*lambda0*x[1]-2.))*(tanh(10*(1. - (lambda0*x[0])-pow(lambda0*x[1],2)/5)) + tanh(2*((lambda0*x[0])+pow(lambda0*x[1],2)/5 + 0.5))) ' \
                  + ' - (x[1]>-3/lambda0 ? 1. : 0.)*(x[1]<=0 ? 1. : 0.)*(lambda0*x[0]>-6 ? 1. : 0.)*ad*0.5*0.5*(1. + tanh(0.5*lambda0*x[1]+2.))*(tanh(10*(1. - (lambda0*x[0])-pow(lambda0*x[1],2)/5)) + tanh(2*((lambda0*x[0])+pow(lambda0*x[1],2)/5 + 0.5))) ' 
        D = seabed
        zeta0 = movingObject
        
        D = Expression(seabed, hd=hd, lambda0=lambda0, element=self.Q.ufl_element())
        self.zeta = Expression(zeta0, ad=ad, c0=c0, t0=self.t0, t=self.t0, hd=hd, lambda0=lambda0, vmax=vmax, element=self.Q.ufl_element())       
        self.h = Expression(seabed + ' + epsilon*(' + movingObject +')', epsilon=epsilon, ad=ad, c0=c0, t0=self.t0, t=self.t0, hd=hd, lambda0=lambda0, vmax=vmax, element=self.Q.ufl_element())
        self.zz_ = interpolate(self.h, self.Q)
        #weak form of the equations
        
        r = 1./dt*inner(U-U_,v)*dx + epsilon*inner(grad(U)*U,v)*dx \
            - div(v)*eta*dx

        r += sigma**2*1./dt*div((D + epsilon*self.zeta)*(U-U_))*div((D + epsilon*self.zeta)*v/2.)*dx \
              - sigma**2*1./dt*div(U-U_)*div((D + epsilon*self.zeta)**2*v/6.)*dx

        r += 1./dt*(eta-eta_)*chi*dx
        r += chi*div((epsilon*eta + D + epsilon*self.zeta)*U)*dx
        #r += chi*(epsilon*eta + D + epsilon*self.zeta)*inner(u,n)*ds

        return r

    def stabilization_parameters(self,U_,eta_,h):
        k1  = 2.
        k2  = 2.
        d1 = k1*(self.dt**(-2) + inner(U_,U_)*h**(-2))**(-0.5)
        d2 = k2*(self.dt**(-2) + eta_*eta_*h**(-2))**(-0.5)

        return d1, d2

    def __str__(self):
          return 'Peregrine'
