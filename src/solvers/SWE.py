__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-28"

from solverbase import *
from numpy import linspace

def f(f0,beta): #Coriolis parameter
    return Expression('f0 + beta*x[1]', f0=f0, beta=beta)

# Class representing the initial conditions
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

# Class for interfacing with the Newton solver
class SWE(NonlinearProblem):
    def __init__(self, problem, W, w, w_k, t, bcs):
        NonlinearProblem.__init__(self)
        h = problem.h
        g = problem.g
        f0 = problem.f0
        beta = problem.beta
        nu = problem.nu

        theta = problem.theta 
        dt = problem.dt

        #define trial and test function
        dw = TrialFunction(W) #direction of the Gateaux derivative
        (v, chi) = TestFunctions(W)

        U, eta = split(w)
        U_k, eta_k = split(w_k)

        # eta_(k+theta)
        eta_mid = (1.0-theta)*eta_k + theta*eta

        # U_(k+theta)
        U_mid = (1.0-theta)*U_k + theta*U

        #weak form of the equations
        L0 = (1/dt)*(eta - eta_k)*chi*dx \
            - h*inner(U_mid,grad(chi))*dx 
        L1 = (1/dt)*inner(U - U_k,v)*dx \
            + f(f0,beta)*(U_mid[0]*v[1] - U_mid[1]*v[0])*dx \
            + g*inner(grad(eta_mid),v)*dx \
            + inner(grad(U_mid)*U_mid,v)*dx \
            + nu*inner(grad(U_mid),grad(v))*dx
        L = L0 + L1

        # Compute directional derivative about w in the direction of dw (Jacobian)
        a = derivative(L, w, dw)

        self.U_k = U_k
        self.eta_k = U_k
        self.L = L
        self.a = a
        self.bcs = bcs
        self.reset_sparsity = True
    def update(self, w_k, bcs, t):
        self.U_k, self.eta_k = split(w_k)
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
        n = int(T/dt)
        t_range = linspace(t, T, n + 1)[1:]

        # Define function spaces
        V = VectorFunctionSpace(mesh, 'CG', 1)
        Q = FunctionSpace(mesh, 'CG', 1)
        W = MixedFunctionSpace([V, Q])

        # Get boundary conditions
        bcs = problem.boundary_conditions(V, Q, t)

        w = Function(W)
        w0 = Function(W)
        U, p = split(w)

        #initial condition
        w0 = InitialConditions(problem, V, Q) 
        w0 = project(w0,W)

        #define the problem
        SWE = SWE(problem, W, w, w0, t, bcs) #build the Shallow Water Equations FE

        #initialize plot
#        viz = plot(w0.split()[1], mesh=mesh, title='Height')

        # Time loop
        self.start_timing()
        for t in t_range:
            #evaluate bcs again (in case they are time-dependent)
            bcs = problem.boundary_conditions(V, Q, t)

            SWE.update(w0, bcs, t) #build the Shallow Water Equations FE

            NewtonSolver(problem.solver).solve(SWE, w.vector()) #solve our problem

            U = w.split()[0]
            eta = w.split()[1]

            # Update
            self.update(problem, t, U, eta)
            #set the solution for the previous time step
            w0.vector()[:] = w.vector() 

            # Plot solution and mesh
#            viz.plot(w0.split()[1]) 
        
        return U, eta
    def __str__(self):
          return 'SWE'

