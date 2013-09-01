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

class Solver(SolverBase):
#    Incremental pressure-correction scheme.

    def __init__(self, options):
        SolverBase.__init__(self, options)

    def solve(self, problem):
        #get problem parameters
        mesh = problem.mesh
        t = 0
        nu = problem.nu #viscosity
        rho = problem.rho #density

        # Define function spaces
        V = VectorFunctionSpace(mesh, 'CG', problem.Pu)
        Q = FunctionSpace(mesh, 'CG', problem.Pp)
        W = MixedFunctionSpace([V, Q])

        # Get boundary conditions
        bcs = problem.boundary_conditions(V, Q, t)

        #define trial and test function
        w = Function(W)
        v, q = TestFunctions(W)
        (U, p) = TrialFunctions(W)

        #weak form of the equations
        L0 = nu*inner(grad(U),grad(v))*dx + 1./rho*inner(grad(p),v)*dx 
        L1 = div(U)*q*dx 
        L = L0 + L1

        # Time loop
        self.start_timing()

        solve(L==0, w, bcs)

        U = w.split()[0]
        p = w.split()[1]

        # Update
        self.update(problem, t, U, p)

        return U, p

    def __str__(self):
          return 'SteadyStokes'
