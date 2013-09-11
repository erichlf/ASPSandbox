__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"

from solverbase import *

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
        W = V * Q

        # Get boundary conditions
        bcs = problem.boundary_conditions(W.sub(0), W.sub(1), t)

        #define trial and test function
        w = Function(W)
        v, q = TestFunctions(W)
        U, p = (as_vector((w[0], w[1])), w[2])

        f = problem.F 
        #weak form of the equations
        F = nu*inner(grad(U),grad(v))*dx + 1./rho*inner(grad(p),v)*dx 
        F -= inner(f(t),v)*dx
        F += div(U)*q*dx 

        # Time loop
        self.start_timing()

        solve(F==0, w, bcs=bcs)

        U, p = w.split()

        # Update
        self.update(problem, t, U, p)

        return U, p

    def __str__(self):
          return 'SteadyStokes'
