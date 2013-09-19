__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"

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
        h = CellSize(mesh) #mesh size

        t = 0
        T = problem.T
        dt = problem.dt
        theta = problem.theta #time stepping method
        nu = problem.nu #viscosity
        rho = problem.rho #density

        # Define function spaces
        V = VectorFunctionSpace(mesh, 'CG', problem.Pu)
        Q = FunctionSpace(mesh, 'CG', problem.Pp)
        W = V * Q

        # Get boundary conditions
        bcs = problem.boundary_conditions(W.sub(0), W.sub(1), t)

        #define trial and test function
        v, q = TestFunctions(W)

        w = Function(W)
        w_ = Function(W)

        #initial condition
        w = InitialConditions(problem, V, Q) 
        w = project(w, W)

        w_ = InitialConditions(problem, V, Q) 
        w_ = project(w_, W)

        U, p = (as_vector((w[0], w[1])), w[2])
        U_, p_ = (as_vector((w_[0], w_[1])), w_[2])

        #U_(k+theta)
        U_theta = (1.0-theta)*U_ + theta*U

        #p_(k+theta)
        p_theta = (1.0-theta)*p_ + theta*p

        f = problem.F
        #weak form of the equations
        F = (1./dt)*inner(U - U_,v)*dx \
            + inner(grad(U_theta)*U_theta,v)*dx \
            + nu*inner(grad(U_theta),grad(v))*dx \
            + 1./rho*inner(grad(p_theta),v)*dx 
        F -= inner(theta*f(t) + (1. - theta)*f(t+theta),v)*dx
        F += div(U_theta)*q*dx 
        if(problem.stabilize):
          # Stabilization parameters
          k1  = 0.5
          k2  = 1.0
          d1 = k1*(dt**(-2) + inner(U_,U_)*h**(-2))**(-0.5) 
          d2 = k2*h 
          #add stabilization
          F += d1*inner(grad(U_theta)*U_theta + grad(p_theta), \
              grad(v)*U_theta + grad(q))*dx + d2*div(U_theta)*div(v)*dx

        U_, p_ = self.timeStepper(problem, t, T, dt, W, w, w_, U_, p_, F) 
        return U_, p_

    def __str__(self):
          return 'NSE'
