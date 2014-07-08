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
        self.Ro = None

    def strong_residual(self,u,U,eta):
        #Parameters
        sigma = problem.sigma
        epsilon = problem.epsilon

        k = problem.k

        zeta_tt = 1./k**2*(self.zeta - 2*self.zeta_ + self.zeta__)
        zeta_t = 1./k*(self.zeta - self.zeta_)

        #strong form for stabilization
        R1 = epsilon*grad(u)*U + grad(eta)
        z1 = -sigma**2/2.*epsilon*self.zeta*grad(zeta_tt)

        R2 = div(u)*(epsilon*eta) + inner(u,grad(epsilon*eta)) \
            + div(U)*epsilon*self.zeta + inner(U,grad(epsilon*self.zeta))
        z2 = zeta_t

        return R1, R2, z1, z2

    def weak_residual(self,W,w,w_,wt,ei_mode=False):
        (U, eta) = (as_vector((w[0], w[1])), w[2])
        (U_, eta_) = (as_vector((w_[0], w_[1])), w_[2])
        (v, chi) = (as_vector((wt[0], wt[1])), wt[2])

        problem = self.problem

        h = CellSize(W.mesh()) #mesh size
        d = 0.1*h**(3./2.) #stabilization parameter

        #set up error indicators
        Z = FunctionSpace(W.mesh(), "DG", 0)
        z = TestFunction(Z)

        alpha = self.alpha #time stepping method

        #Parameters
        sigma = problem.sigma
        epsilon = problem.epsilon

        t0 = problem.t0 #initial time
        self.t0 = t0
        T = problem.T #Final time
        k = problem.k #time step

        D, self.zeta, self.zeta_, self.zeta__, self.bottom, self.H, self.H_ \
            = self.seabed(problem,self.Q,t0,epsilon)

        zeta_tt = 1./k**2*(self.zeta - 2*self.zeta_ + self.zeta__)
        zeta_t = 1./k*(self.zeta - self.zeta_)

        #Time stepping method
        U_alpha = (1. - alpha)*U_ + alpha*U
        eta_alpha = (1. - alpha)*eta_ + alpha*eta
        zeta_alpha = (1. - alpha)*self.zeta_ + alpha*self.zeta

        t = t0 + k
        #forcing and mass source/sink
        F1_alpha = alpha*problem.F1(t) + (1 - alpha)*problem.F1(t0)
        F2_alpha = alpha*problem.F2(t) + (1 - alpha)*problem.F2(t0)


        if(not self.options["stabilize"] or ei_mode):
          d = 0
        if(not ei_mode):
          z = 1.

        #weak form of the equations
        r = z*(1./k*inner(U-U_,v) + epsilon*inner(grad(U_alpha)*U_alpha,v) \
            - div(v)*eta_alpha)*dx

        r += z*(sigma**2*1./k*div((D + epsilon*zeta_alpha)*(U-U_))*div((D + epsilon*zeta_alpha)*v/2.) \
              - sigma**2*1./k*div(U-U_)*div((D + epsilon*zeta_alpha)**2*v/6.))*dx
        r += z*sigma**2*zeta_tt*div((D + epsilon*zeta_alpha)*v/2.)*dx

        r += z*(1./k*(eta-eta_)*chi + zeta_t*chi)*dx
        r -= z*inner(U_alpha,grad(chi))*(epsilon*eta_alpha + D + epsilon*zeta_alpha)*dx

        r -= z*(inner(F1_alpha,v) + F2_alpha*chi)*dx

        r += z*d*(inner(grad(U_alpha),grad(v)) + inner(grad(eta_alpha),grad(chi)))*dx

        return r

    def stabilization_parameters(self,U_,eta_,h):
        K1  = 2.
        K2  = 2.
        d1 = K1*(self.k**(-2) + inner(U_,U_)*h**(-2))**(-0.5)
        d2 = K2*(self.k**(-2) + eta_*eta_*h**(-2))**(-0.5)

        return d1, d2

    def functional(self,mesh,u,eta):

      M = u[0]*dx # Mean of the x-velocity in the whole domain

      return M

    def seabed(self,problem,Q,t0,epsilon):
        D = Expression(problem.D, element=Q.ufl_element())
        zeta = Expression(problem.zeta0, t0=t0, t=t0,element=Q.ufl_element())
        zeta_ = Expression(problem.zeta0, t0=t0, t=t0, element=Q.ufl_element())
        zeta__ = Expression(problem.zeta0, t0=t0, t=t0, element=Q.ufl_element())
        bottom = Expression(problem.D + ' + epsilon*(' + problem.zeta0 +')', \
                epsilon=epsilon, t0=t0, t=t0, element=Q.ufl_element())
        H = interpolate(bottom, Q)
        H_ = interpolate(bottom, Q)

        return D, zeta, zeta_, zeta__, bottom, H, H_

    def wave_object(self, Q, t, k):
        self.zeta.t = t
        self.zeta_.t = max(t - k, self.t0)
        self.zeta__.t = max(t - 2*k, self.t0)
        self.H_.assign(self.H)
        self.bottom.t = t
        self.H = interpolate(self.bottom, Q)

    def __str__(self):
          return 'Peregrine'
