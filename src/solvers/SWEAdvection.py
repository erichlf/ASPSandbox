__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"

from solverbase import *

def f(f0,beta): #Coriolis parameter
    return Expression('f0 + beta*x[1]', f0=f0, beta=beta)

class InitialConditions(Expression):
    def __init__(self, problem, V, Q, S):
        self.U0, self.eta0 = problem.initial_conditions(V, Q)
        self.U0 = project(self.U0,V)
        self.eta0 = project(self.eta0,Q)
        self.s0 = Expression('A*exp(-(pow(x[0]-0.75,2)+pow(x[1]-0.75,2))/(2*S*S))', A=1, S=0.05)
        self.s0 = project(self.s0,S)

    def eval(self, value, x):
        value[:2] = self.U0(x)
        value[2] = self.eta0(x)
        value[3] = self.s0(x)

    def value_shape(self):
        return (4,)

class Solver(SolverBase):
#    Incremental pressure-correction scheme.

    def __init__(self, options):
        SolverBase.__init__(self, options)

    def solve(self, problem):
        #get problem mesh 
        mesh = problem.mesh

        t = 0 #initial time
        T = problem.T #final time
        dt = problem.dt #time step
        theta = problem.theta #time stepping method

        #problem parameters
        h = problem.h #fluid depth
        g = problem.g #gravity
        f0 = problem.f0 #reference Coriolis force
        beta = problem.beta #beta plane parameter
        nu = problem.nu #viscosity
        rho = problem.rho #density
        K = 1E-3 #diffusion coefficient

        # Define function spaces
        V = VectorFunctionSpace(mesh, 'CG', problem.Pu)
        Q = FunctionSpace(mesh, 'CG', problem.Pp)
        S = FunctionSpace(mesh, 'CG', problem.Pp)
        W = MixedFunctionSpace([V, Q, S])

        #define trial and test function
        v, chi, r = TestFunctions(W)

        w = Function(W)
        w_ = Function(W)

        #initial condition
        w = InitialConditions(problem, V, Q, S) 
        w = project(w, W)

        w_ = InitialConditions(problem, V, Q, S) 
        w_ = project(w_, W)

        U, eta, s = (as_vector((w[0], w[1])), w[2], w[3])
        U_, eta_, s_ = (as_vector((w_[0], w_[1])), w_[2], w_[3])

        # Get boundary conditions
        bcs = problem.boundary_conditions(W.sub(0), W.sub(1), t)

        #U_(k+theta)
        U_theta = (1.0-theta)*U_ + theta*U

        #p_(k+theta)
        eta_theta = (1.0-theta)*eta_ + theta*eta

        #s_(k+theta)
        s_theta = (1.0-theta)*s_ + theta*s

        #weak form of the equations
        F = (1./dt)*(eta - eta_)*chi*dx \
            + h*div(U_theta)*chi*dx 
        F += (1./dt)*inner(U - U_,v)*dx \
            + f(f0,beta)*(U_theta[0]*v[1] - U_theta[1]*v[0])*dx \
            - g*eta_theta*div(v)*dx \
            + inner(grad(U_theta)*U_theta,v)*dx \
            + nu*inner(grad(U_theta),grad(v))*dx
        F += (1./dt)*(s - s_)*r*dx \
            + inner(U_theta,grad(s_theta))*r*dx \
            + K*inner(grad(s_theta),grad(r))*dx

        # Time loop
        self.start_timing()

        #plot and save initial condition
        self.update(problem, t, w_.split()[0], w_.split()[2]) 

        while t<T:
            t += dt

            #evaluate bcs again (in case they are time-dependent)
            bcs = problem.boundary_conditions(W.sub(0), W.sub(1), t)

            solve(F==0, w, bcs=bcs)

            w_.vector()[:] = w.vector()

            U_ = w_.split()[0] 
            eta_ = w_.split()[1]
            s_ = w_.split()[2]

            # Update
            self.update(problem, t, U_, s_)

        return U_, s_

    def __str__(self):
          return 'SWEAdvection'
