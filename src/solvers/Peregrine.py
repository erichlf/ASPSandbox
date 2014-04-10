__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"

from solverbase import *

class Solver(SolverBase):
#    Incremental pressure-correction scheme.

    def __init__(self, options):
        SolverBase.__init__(self, options)

        #parameters
        self.sigma = 1./options['Re']
        self.eps = options['Theta']

    def wave_object(self,W,t):
        """
            Here we create a moving object under water to create a wave.
            Since we are using a mixed space we have to project by solving
            a system. Then we can extract zeta from the solution to the
            system.
        """

        #zeta = Expression('A/sqrt(2*S*pi)*exp(-pow(x[0] - (t - 0.5),2)/(2*S))',
        #L is the length of half the domain
        #P is the amplitude
        #U is the speed of the moving object
        Lx = (self.x1 - self.x0)
        Ly = (self.y1 - self.y0)

        zeta0 = 'a0*exp(-pow(x[0]-xh-vh*t,2)/(bh*bh))'
        zeta = Expression(zeta0, a0=0.1, xh=0.3, vh=1, bh=0.02,t=t)

        zeta = self.Q_project(zeta,W)

        return zeta

    def strong_residual(self,u,U,eta,zeta,zeta_t,zeta_tt):
        eps = self.eps

        R1 = (1. + zeta + eps*eta)*div(U) \
            + inner(grad(1 + zeta + eps*eta),U)
        R1 += zeta_t
        R2 = grad(eta)
        R2 -= (1. + zeta)/2.*grad(zeta_tt)

        return R1, R2

    def weak_residual(self,U,U_,eta,eta_,zeta,zeta_,zeta__,v,chi):
        dt = self.dt #time step
        alpha = self.alpha #time stepping method

        #parameters
        eps = self.eps
        sigma = self.sigma

        #U_(k+alpha)
        U_alpha = (1.0-alpha)*U_ + alpha*U

        #p_(k+alpha)
        eta_alpha = (1.0-alpha)*eta_ + alpha*eta

        zeta_t = 1./dt*(zeta - zeta_)
        zeta_tt = 1./dt**2*(zeta - 2*zeta_ + zeta__)

        #weak form of the equations
        r = (1./dt)*(eta - eta_)*chi*dx \
            + (1. + zeta + eps*eta_alpha)*div(U_alpha)*chi*dx \
            + inner(grad(1 + zeta + eps*eta_alpha),U_alpha)*chi*dx
        r += zeta_t*chi*dx
        r += (1./dt)*inner(U - U_,v)*dx \
            + inner(grad(eta_alpha),v)*dx \
            + eps*inner(grad(U_alpha)*U_alpha,v)*dx \
            + 1./dt*sigma**2*(1. + zeta)**2/3.*inner(grad(U - U_),grad(v))*dx
        r -= (1. + zeta)/2.*inner(grad(zeta_tt),v)*dx

        return r

    def stabilization_parameters(self,U_,eta_,h):
        k1  = self.eps/2
        k2  = 1/2
        d1 = k1*(self.dt**(-2) + inner(U_,U_)*h**(-1))**(-0.5)
        d2 = k2*(self.dt**(-2) + eta_*eta_*h**(-1))**(-0.5)

        return d1, d2

    def solve(self, problem):
        #get problem mesh
        mesh = problem.mesh
        h = CellSize(mesh) #mesh size
        self.x0 = problem.mesh.coordinates()[0][0]
        self.x1 = problem.mesh.coordinates()[-1][0]
        self.y0 = problem.mesh.coordinates()[0][1]
        self.y1 = problem.mesh.coordinates()[-1][1]

        t = self.t0 #initial time
        T = problem.T #final time
        dt = self.dt #time step
        alpha = self.alpha #time stepping method

        Pu = self.Pu
        Pp = self.Pp

        # Define function spaces
        V = VectorFunctionSpace(mesh, 'CG', Pu)
        Q = FunctionSpace(mesh, 'CG', Pp)
        W = MixedFunctionSpace([V, Q])

        # Get boundary conditions
        bcs = problem.boundary_conditions(W.sub(0), W.sub(1), t)

        #define trial and test function
        v, chi = TestFunctions(W)

        #initialize w as a function in W
        w = Function(W)
        w_ = Function(W)

        #initialize zetas as functions in Q
        zeta = w.split()[1]
        zeta_ = w.split()[1]
        zeta__ = w.split()[1]

        #initial condition
        w = self.InitialConditions(problem, W)
        w_ = self.InitialConditions(problem, W)

        #assign zeta as the wave object
        zeta = self.wave_object(W,t)
        zeta_.assign(zeta)
        zeta__.assign(zeta)

        #initialize U and eta as functions in V and Q, respectively
        U, eta = (as_vector((w[0], w[1])), w[2])
        U_, eta_ = (as_vector((w_[0], w_[1])), w_[2])

        #U_(k+alpha)
        U_alpha = (1.0-alpha)*U_ + alpha*U

        #p_(k+alpha)
        eta_alpha = (1.0-alpha)*eta_ + alpha*eta

        zeta_t = 1./dt*(zeta - zeta_)
        zeta_tt = 1./dt**2*(zeta - 2*zeta_ + zeta__)

        F = self.weak_residual(U, U_, eta, eta_, zeta, zeta_, zeta__, v, chi)

        U_, eta_ = self.timeStepper(problem, t, T, self.dt, W, w, w_, U_, eta_, zeta, zeta_, zeta__, F)
        return U_, eta_

    def timeStepper(self, problem, t, T, dt, W, w, w_, U_, eta_, zeta, zeta_, zeta__, F):
        # Time loop
        self.start_timing()

        #plot and save initial condition
        self.update(problem, t, zeta, w_.split()[1])
        while t<T:
            t += dt
            #update the wave generator
            zeta__.assign(zeta_)
            zeta_.assign(zeta)
            zeta = self.wave_object(W,t)

            #evaluate bcs again (in case they are time-dependent)
            bcs = problem.boundary_conditions(W.sub(0), W.sub(1), t)

            solve(F==0, w, bcs=bcs)

            w_.vector()[:] = w.vector()

            U_ = w_.split()[0]
            eta_ = w_.split()[1]

            # Update
            self.update(problem, t, zeta, eta_)

        return U_, eta_

    def __str__(self):
          return 'Peregrine'
