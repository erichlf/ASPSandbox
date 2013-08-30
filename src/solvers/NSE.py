__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"

from solverbase import *
from numpy import linspace

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
    def __init__(self, problem, W, w, w_k, t, bcs):
        NonlinearProblem.__init__(self)

        theta = problem.theta #time stepping method
        dt = problem.dt #time step
        nu = problem.nu #viscosity
        rho = problem.rho #density

        #define trial and test function
        dw = TrialFunction(W) #direction of the Gateaux derivative
        v, q = TestFunctions(W)

        U, p = split(w)
        U_k, p_k = split(w_k)

        #U_(k+theta)
        U_mid = (1.0-theta)*U_k + theta*U

        #p_(k+theta)
        p_mid = (1.0-theta)*p_k + theta*p

        #weak form of the equations
        L0 = (1./dt)*inner(U - U_k,v)*dx \
            + inner(grad(U_mid)*U_mid,v)*dx \
            + nu*inner(grad(U_mid),grad(v))*dx \
            - 1./rho*p_mid*div(v)*dx 
        L1 = div(U_mid)*q*dx 
        L = L0 + L1

        # Compute directional derivative about w in the direction of dw (Jacobian)
        a = derivative(L, w, dw)

        self.t = t
        self.U_k = U_k
        self.p_k = p_k
        self.L = L
        self.a = a
        self.bcs = bcs
        self.reset_sparsity = True
    def update(self, w_k, bcs, t):
        self.U_k, self.p_k = split(w_k)
        self.bcs = bcs
        self.t = t
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

        # Define function spaces
        V = VectorFunctionSpace(mesh, 'CG', problem.Pu)
        Q = FunctionSpace(mesh, 'CG', problem.Pp)
        W = MixedFunctionSpace([V, Q])

        # Get boundary conditions
        bcs = problem.boundary_conditions(V, Q, t)

        w = Function(W)
        w0 = Function(W)
        U, p = split(w)

        #initial condition
        w0 = InitialConditions(problem, V, Q) 
        w0 = project(w0,W)

        # Time loop
        self.start_timing()

        #plot and save initial condition
        self.update(problem, t, w0.split()[0], w0.split()[1]) 

        #create problem 
        NS = NSE(problem, W, w, w0, t, bcs) #build problem 

        while t<T:
            t += dt

            #evaluate bcs again (in case they are time-dependent)
            bcs = problem.boundary_conditions(V, Q, t)

            NS.update(w0, bcs, t)
            NewtonSolver(problem.solver).solve(NS, w.vector())

            U = w.split()[0]
            p = w.split()[1]

            # Update
            self.update(problem, t, U, p)
            #set the solution for the previous time step
            w0.vector()[:] = w.vector() 
        
        return U, p

    def __str__(self):
          return 'NSE'
