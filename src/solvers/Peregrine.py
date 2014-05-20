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
        bh = 0.7 #width of the moving object

        #Scaled Parameters
        self.dt = self.dt*c0/lambda0 #time step
        dt = self.dt
        self.T = self.T*c0/lambda0 #Final time
        hd = hd/h0 #depth
        ad = ad/a0 #height of the moving object

        #Definition of the wave_object
        vmax = ((hd*h0+ad*a0)*g)**(0.5) #Max Speed of the moving object [m.s^(-1)]
        seabed = 'hd - 0.5/10.*(x[1]>4./lambda0 ? 1. : 0.)*(lambda0*x[1]-4.) + 0.5/10.*(x[1]<(-4./lambda0) ? 1. : 0.)*(lambda0*x[1]+4.)'
        #traj = 'vmax*lambda0/c0*t*exp(-0.001/pow(lambda0/c0*t,2))'
        traj = 'vmax*lambda0/c0*t*exp(-4./(lambda0/c0*t+0.05))'
        movingObject = ' - (x[1]<3/lambda0 ? 1. : 0.)*(x[1]>0 ? 1. : 0.)*(lambda0*x[0]-'+traj+'>-6 ? 1. : 0.)*epsilon*ad*0.5*0.5*(1. - tanh(0.5*lambda0*x[1]-2.))*(tanh(10*(1. - (lambda0*x[0] - ' + traj + ')-pow(lambda0*x[1],2)/5)) + tanh(2*((lambda0*x[0] - ' + traj + ')+pow(lambda0*x[1],2)/5 + 0.5))) ' \
                  + ' - (x[1]>-3/lambda0 ? 1. : 0.)*(x[1]<=0 ? 1. : 0.)*(lambda0*x[0]-'+traj+'>-6 ? 1. : 0.)*epsilon*ad*0.5*0.5*(1. + tanh(0.5*lambda0*x[1]+2.))*(tanh(10*(1. - (lambda0*x[0] - ' + traj + ')-pow(lambda0*x[1],2)/5)) + tanh(2*((lambda0*x[0] - ' + traj + ')+pow(lambda0*x[1],2)/5 + 0.5))) ' 
        zeta0 = seabed + movingObject
        
        self.zeta = Expression(zeta0, ad=ad, c0=c0, t0=self.t0, t=self.t0, bh=bh, hd=hd, epsilon=epsilon, lambda0=lambda0, vmax=vmax, element=self.Q.ufl_element())
        self.zeta_ = Expression(zeta0, ad=ad, c0=c0, t0=self.t0, t=self.t0, bh=bh, hd=hd, epsilon=epsilon, lambda0=lambda0, vmax=vmax, element=self.Q.ufl_element())
        self.zeta__ = Expression(zeta0, ad=ad, c0=c0, t0=self.t0, t=self.t0, bh=bh, hd=hd, epsilon=epsilon, lambda0=lambda0, vmax=vmax, element=self.Q.ufl_element())

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
        self.zz_.assign(self.z)
        self.z = interpolate(self.zeta, self.Q)

    def __str__(self):
          return 'Peregrine'
