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

class NSE(NonlinearProblem):
    def __init__(self, a, L, bcs):
        NonlinearProblem.__init__(self)

        self.L = L
        self.a = a
        self.bcs = bcs
        self.reset_sparsity = True
    def F(self, b, x):
        assemble(self.L, tensor=b, bcs=self.bcs)
    def J(self, A, x):
        assemble(self.a, tensor=A, reset_sparsity=self.reset_sparsity,
            bcs=self.bcs)
        #self.reset_sparsity = False

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
        dw = TrialFunction(W) #direction of the Gateaux derivative
        v, q = TestFunctions(W)

        w = Function(W)
        w_ = Function(W)

        #initial condition
        w = InitialConditions(problem, V, Q) 
        w = project(w, W)

        w_ = InitialConditions(problem, V, Q) 
        w_ = project(w_, W)

        U, p = split(w)
        U_, p_ = split(w_)

        #U_(k+theta)
        U_mid = (1.0-theta)*U_ + theta*U

        #p_(k+theta)
        p_mid = (1.0-theta)*p_ + theta*p

        #weak form of the equations
        L0 = (1./dt)*inner(U - U_,v)*dx \
            + inner(grad(U_mid)*U_mid,v)*dx \
            + nu*inner(grad(U_mid),grad(v))*dx \
            + 1./rho*inner(grad(p_mid),v)*dx 
        L1 = div(U_mid)*q*dx 
        L = L0 + L1

        # Compute directional derivative about w in the direction of dw (Jacobian)
        a = derivative(L, w, dw)

        # Time loop
        self.start_timing()

        #plot and save initial condition
        self.update(problem, t, w_.split()[0], w_.split()[1]) 

        #create problem 
        NS = NSE(a, L, bcs) #build problem 

        while t<T:
            t += dt

            w_.vector()[:] = w.vector()

            #evaluate bcs again (in case they are time-dependent)
            bcs = problem.boundary_conditions(V, Q, t)
            NS.bcs = bcs

            NewtonSolver(problem.solver).solve(NS, w.vector())

            # Update
            self.update(problem, t, w.split()[0], w.split()[1])
            #set the solution for the previous time step
        
        return U_, p_

    def __str__(self):
          return 'NSE'
