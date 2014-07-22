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

    #define functions spaces
    def function_space(self, mesh):
        V = VectorFunctionSpace(mesh, 'CG', self.Pu)
        Q = FunctionSpace(mesh, 'CG', self.Pp)
        self.Q = Q #Bad hack for being able to project into Q space
        W = MixedFunctionSpace([V, Q])

        self.Zeta = Function(Q, name='zeta')
        self.Zeta_ = Function(Q, name='zeta_')
        self.Zeta__ = Function(Q, name='zeta__')
        self.zeta0 = Function(Q, name='InitialShape')

        return W

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

    def weak_residual(self,problem, W, w, w_, wt, ei_mode=False):
        (U, eta) = (as_vector((w[0], w[1])), w[2])
        (U_, eta_) = (as_vector((w_[0], w_[1])), w_[2])
        (v, chi) = (as_vector((wt[0], wt[1])), wt[2])

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

        #set up our wave object
        self.wave_object(problem, self.Q, t0, k)

        #save the initial wave object
        self.zeta0.assign(self.Zeta)

        zeta_tt = 1./k**2*(self.Zeta - 2*self.Zeta_ + self.Zeta__)
        zeta_t = 1./k*(self.Zeta - self.Zeta_)

        #Time stepping method
        U_alpha = (1. - alpha)*U_ + alpha*U
        eta_alpha = (1. - alpha)*eta_ + alpha*eta
        zeta_alpha = (1. - alpha)*self.Zeta_ + alpha*self.Zeta
        H_alpha = (1. - alpha)*self.H_ + alpha*self.H

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

        r += z*(sigma**2*1./k*div(H_alpha*(U-U_))*div(H_alpha*v/2.) \
              - sigma**2*1./k*div(U-U_)*div(H_alpha**2*v/6.))*dx
        r += z*sigma**2*zeta_tt*div(H_alpha*v/2.)*dx

        r += z*(1./k*(eta-eta_)*chi + zeta_t*chi)*dx
        r -= z*inner(U_alpha,grad(chi))*(epsilon*eta_alpha + H_alpha)*dx

        r -= z*(inner(F1_alpha,v) + F2_alpha*chi)*dx

        r += z*d*(inner(grad(U_alpha),grad(v)) + inner(grad(eta_alpha),grad(chi)))*dx

        return r

    def stabilization_parameters(self,U_,eta_,h):
        K1  = 2.
        K2  = 2.
        d1 = K1*(self.k**(-2) + inner(U_,U_)*h**(-2))**(-0.5)
        d2 = K2*(self.k**(-2) + eta_*eta_*h**(-2))**(-0.5)

        return d1, d2

    #Optimization Function
    def Optimize(self, problem, w):
        eta = w.split()[1]
        #bounds on object
        lb = project(Expression('-0.5'), self.Q, name='LowerBound')
        ub = project(Expression('0.0'), self.Q, name='UpperBound')

        #Functionnal to be minimized: L2 norme over a subdomain
        J = Functional(-inner(eta, eta)*dx*dt[FINISH_TIME] + self.Zeta*self.Zeta*dx)

        #initial shape
        m = ScalarParameters(problem.params)
        Jhat = ReducedFunctional(J, m) #Reduced Functional
        shape_opt = minimize(Jhat, bounds=(lb, ub))

        return shape_opt

    def functional(self,mesh,w):

        (u, eta) = (as_vector((w[0], w[1])), w[2])

        M = u[0]*dx # Mean of the x-velocity in the whole domain

        return M

    def seabed(self, problem, Q, t, k, epsilon):
        D = problem.D

        problem.zeta0.t = t
        zeta = project(problem.zeta0, Q)
        problem.zeta0.t = max(t - k, self.t0)
        zeta_ = project(problem.zeta0, Q)
        problem.zeta0.t = max(t - 2*k, self.t0)
        zeta__ = project(problem.zeta0, Q)

        H = project(D + epsilon*zeta, Q, name='Bathymetry')
        H_ = project(D + epsilon*zeta_, Q, name='PreviousBathymetry')

        return D, zeta, zeta_, zeta__, H, H_

    def wave_object(self, problem, Q, t, k):
        D, zeta, zeta_, zeta__, self.H, self.H_ \
            = self.seabed(problem, Q, t, k, problem.epsilon)

        self.Zeta.assign(zeta)
        self.Zeta_.assign(zeta_)
        self.Zeta__.assign(zeta__)

    def Plot(self, problem, W, w):
        u = w.split()[0]
        eta = w.split()[1]

        # Plot velocity and height and wave object
        if self.vizU is None:
            self.vizU = plot(u, title='Velocity', rescale=True)
            self.vizP = plot(eta, title='Height', rescale=True)
            self.vizZ = plot(self.H_, title='Wave Object', rescale=True)
        else :
            self.vizU.plot(u)
            self.vizP.plot(eta)
            self.vizZ.plot(self.H_)

    def __str__(self):
          return 'Peregrine'
