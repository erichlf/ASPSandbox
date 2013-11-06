__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"

from solverbase import *

class InitialConditions(Expression):
    def __init__(self, problem, V, Q):
        self.U0, self.eta0 = problem.initial_conditions(V, Q)
        self.U0 = project(self.U0,V)
        self.eta0 = project(self.eta0,Q)

    def eval(self, value, x):
        value[:2] = self.U0(x)
        value[2] = self.eta0(x)

    def value_shape(self):
        return (3,)

def wave_object(Q,t):
    zeta = Expression('A/sqrt(2*S*pi)*exp(-pow(x[0] - (t - 0.5),2)/(2*S))', A=200., S=0.02, t=t)
    zeta = project(zeta,Q)

    return zeta

class Solver(SolverBase):
#    Incremental pressure-correction scheme.

    def __init__(self, options):
        SolverBase.__init__(self, options)

    def solve(self, problem):
        #get problem mesh
        mesh = problem.mesh
        h = CellSize(mesh) #mesh size

        t = 0 #initial time
        T = problem.T #final time
        dt = problem.dt #time step
        theta = problem.theta #time stepping method

        #problem parameters
        H = problem.h #fluid depth
        g = problem.g #gravity
        f0 = problem.f0 #reference Coriolis force
        beta = problem.beta #beta plane parameter
        nu = problem.nu #viscosity
        rho = problem.rho #density

        eps = nu
        sigma = nu

        # Define function spaces
        V = VectorFunctionSpace(mesh, 'CG', problem.Pu)
        Q = FunctionSpace(mesh, 'CG', problem.Pp)
        W = V * Q

        # Get boundary conditions
        bcs = problem.boundary_conditions(W.sub(0), W.sub(1), t)

        #define trial and test function
        v, chi = TestFunctions(W)

        w = Function(W)
        w_ = Function(W)

        #initial condition
        w = InitialConditions(problem, V, Q)#V, Q)
        #w = project(w, W)

        w_ = InitialConditions(problem, W.sub(0), W.sub(1))#V, Q)
        #w_ = project(w_, W)

        zeta = wave_object(W.sub(1),t)
        zeta_ = project(zeta,W.sub(1))
        zeta__ = project(zeta,W.sub(1))

        U, eta = (as_vector((w[0], w[1])), w[2])
        U_, eta_ = (as_vector((w_[0], w_[1])), w_[2])

        #U_(k+theta)
        U_theta = (1.0-theta)*U_ + theta*U

        #p_(k+theta)
        eta_theta = (1.0-theta)*eta_ + theta*eta

        #weak form of the equations
        F = (1./dt)*(eta - eta_)*chi*dx \
            + (H + eps*eta_theta)*div(U_theta)*chi*dx \
            + inner(grad(H + eps*eta_theta),U_theta)*chi*dx
        F += 1./dt*(zeta - zeta_)*chi*dx
        F += (1./dt)*inner(U - U_,v)*dx \
            - eta_theta*div(v)*dx \
            + eps*inner(grad(U_theta)*U_theta,v)*dx \
            + 1./dt*sigma**2*H**2/3.*inner(grad(U - U_),grad(v))*dx
        F -= H/2.*inner(grad(1./dt**2*(zeta - 2*zeta_ + zeta__)),v)*dx

        if(problem.stabilize):
          # Stabilization parameters
          k1  = 0.5
          k2  = 0.5
          d1 = k1*(dt**(-2) + inner(U_,U_)*h**(-2))**(-0.5)
          d2 = k2*(dt**(-2) + eta_*eta_*h**(-2))**(-0.5)

          #add stabilization
          F += d1*(inner(grad(H+eps*eta_theta),U_theta) \
                  + (H + eps*eta_theta)*div(U_theta) + 1./dt*(zeta - zeta_)) \
              *(inner(grad(H+eps*eta_theta),v) + (H + eps*eta_theta)*div(v))*dx
          F += d2*inner(grad(eta_theta) + eps*grad(U_theta)*U_theta \
                  + H/2.*grad(1./dt**2*(zeta - 2*zeta_ + zeta__)), \
              grad(chi) + eps*grad(v)*U_theta)*dx

        # Time loop
        self.start_timing()

        #plot and save initial condition
        self.update(problem, t, w_.split()[0], w_.split()[1])
        #zt = project(1./dt*(zeta - zeta_),Q)
        #viz = plot(zt, rescale=True)

        while t<T:
            t += dt
            #update the wave generator
            zeta__.assign(zeta_)
            zeta_.assign(zeta)
            zeta = wave_object(W.sub(1),t)
            #zt = project(1./dt*(zeta - zeta_),Q)
            #viz.update(zt)

            #evaluate bcs again (in case they are time-dependent)
            bcs = problem.boundary_conditions(W.sub(0), W.sub(1), t)

            solve(F==0, w, bcs=bcs)

            w_.vector()[:] = w.vector()

            U_ = w_.split()[0]
            eta_ = w_.split()[1]

            # Update
            self.update(problem, t, U_, eta_)

        return U_, eta_

    def __str__(self):
          return 'MovingSWE'
