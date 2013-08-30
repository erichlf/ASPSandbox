__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"

from solverbase import *

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
class LinearSWE(NonlinearProblem):
    def __init__(self, problem, W, w, w_k, t, bcs):
        NonlinearProblem.__init__(self)
        h = problem.h
        g = problem.g
        f0 = problem.f0
        beta = problem.beta

        theta = problem.theta 
        dt = problem.dt

        #define trial and test function
        (v, chi) = TestFunctions(W)

        U, eta = split(w)
        U_k, eta_k = split(w_k)

        # eta_(k+theta)
        eta_mid = (1.0-theta)*eta_k + theta*eta

        # U_(k+theta)
        U_mid = (1.0-theta)*U_k + theta*U

        #weak form of the equations
        L0 = (1./dt)*(eta - eta_k)*chi*dx \
            + h*div(U_mid)*chi*dx 
        L1 = (1./dt)*inner(U - U_k,v)*dx \
            + f(f0,beta)*(U_mid[0]*v[1] - U_mid[1]*v[0])*dx \
            + g*inner(grad(eta_mid),v)*dx 
        L = L0 + L1

        # Compute directional derivative about w in the direction of dw (Jacobian)

        self.U_k = U_k
        self.eta_k = U_k
        self.L = L
        self.bcs = bcs
        self.reset_sparsity = True
    def update(self, w_k, bcs, t):
        self.U_k, self.eta_k = split(w_k)
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

        #define the problem
        SWE = LinearSWE(problem, W, w, w0, t, bcs) #build the Shallow Water Equations FE

        while t<T:
            t += dt

            #evaluate bcs again (in case they are time-dependent)
            bcs = problem.boundary_conditions(V, Q, t)

            SWE.update(w0, bcs, t) #build the Shallow Water Equations FE

            solve(SWE.L==0, w, bcs) #solve our problem

            U = w.split()[0]
            eta = w.split()[1]

            # Update
            self.update(problem, t, U, eta)
            #set the solution for the previous time step
            w0.vector()[:] = w.vector() 

        return U, eta

    def __str__(self):
          return 'LinearSWE'

