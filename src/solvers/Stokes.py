__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"

from solverbase import *

class InitialConditions(Expression):
    def __init__(self, problem, V, Q):
        self.U0, self.p0 = problem.initial_conditions(V, Q)
        self.U0 = project(self.U0,V)
        self.p0 = project(self.p0,Q)

    def eval(self, value, x):
        value[:2] = self.U0(x)
        value[2] = self.p0(x)

    def value_shape(self):
        return (3,)

class Stokes(NonlinearProblem):
    def __init__(self, problem, W, w, w_k, t, bcs):
        NonlinearProblem.__init__(self)

        theta = problem.theta #time stepping method
        dt = problem.dt #time step
        nu = problem.nu #viscosity
        rho = problem.rho #density

        #define trial and test function
        v, q = TestFunctions(W)

        U, p = split(w)
        U_k, p_k = split(w_k)

        #U_(k+theta)
        U_mid = (1.0-theta)*U_k + theta*U

        #p_(k+theta)
        p_mid = (1.0-theta)*p_k + theta*p

        #weak form of the equations
        L0 = (1./dt)*inner(U - U_k,v)*dx \
            + nu*inner(grad(U_mid),grad(v))*dx \
            + 1./rho*inner(grad(p_mid),v)*dx 
        L1 = inner(div(U_mid),q)*dx #negative makes things symmetric 
        L = L0 + L1

        self.t = t
        self.U_k = U_k
        self.p_k = p_k
        self.L = L
        self.bcs = bcs
        self.reset_sparsity = True
    def update(self, w_k, bcs, t):
        self.U_k, self.p_k = split(w_k)
        self.bcs = bcs
        self.t = t

class Solver(SolverBase):
#    Incremental pressure-correction scheme.

    def __init__(self, options):
        SolverBase.__init__(self, options)

    def solve(self, problem):
        #get problem parameters
        mesh = problem.mesh
        t = 0
        T = problem.T
        dt = problem.dt
        theta = problem.theta #time stepping method
        nu = problem.nu #viscosity
        rho = problem.rho #density

        # Define function spaces
        V = VectorFunctionSpace(mesh, 'CG', problem.Pu)
        Q = FunctionSpace(mesh, 'CG', problem.Pp)
        W = MixedFunctionSpace([V, Q])

        # Get boundary conditions
        bcs = problem.boundary_conditions(V, Q, t)

        #define trial and test function
        v, q = TestFunctions(W)
        U, p = TrialFunctions(W)

        w = Function(W)
        w_ = Function(W)

        #initial condition
        w = InitialConditions(problem, V, Q) 
        w = project(w,W)

        w_ = InitialConditions(problem, V, Q) 
        w_ = project(w_,W)

        U_, p_ = split(w_)

        #U_(k+theta)
        U_mid = (1.0-theta)*U_ + theta*U

        #p_(k+theta)
        p_mid = (1.0-theta)*p_ + theta*p

        #weak form of the equations
        F0 = (1./dt)*inner(U - U_,v)*dx \
            + nu*inner(grad(U_mid),grad(v))*dx \
            + 1./rho*inner(grad(p_mid),v)*dx 
        F1 = div(U_mid)*q*dx 
        F = F0 + F1
        a = lhs(F)
        L = rhs(F)

        # Time loop
        self.start_timing()

        #plot and save initial condition
        self.update(problem, t, w_.split()[0], w_.split()[1]) 

        while t<T:
            t += dt

            w_.vector()[:] = w.vector()

            U_ = w_.split()[0]
            p_ = w_.split()[1]

            #evaluate bcs again (in case they are time-dependent)
            bcs = problem.boundary_conditions(V, Q, t)

            A = assemble(a)
            b = assemble(L)
            if isinstance(bcs, list):
                [bc.apply(A, b) for bc in bcs]
            else :
                bcs.apply(A, b)
            solve(A, w.vector(), b)

            # Update
            self.update(problem, t, w.split()[0], w.split()[1])
        
        return U_, p_

    def __str__(self):
          return 'Stokes'
