__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"

from solverbase import *

class InitialConditions(Expression):
    def __init__(self, problem, V, Q, R):
        A = 1
        sig = 0.05
        self.U0, self.p0 = problem.initial_conditions(V, Q)
        self.U0 = project(self.U0,V)
        self.p0 = project(self.p0,Q)
        self.rho0 = Expression('A*exp(-(pow(x[0]-0.75,2)+pow(x[1]-0.75,2))/(2*sig*sig))', A=A, sig=sig)
        self.rho0 = project(self.rho0,R)

    def eval(self, value, x):
        value[:2] = self.U0(x)
        value[2] = self.p0(x)
        value[3] = self.rho0(x)

    def value_shape(self):
        return (4,)


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

        # Define function spaces
        V = VectorFunctionSpace(mesh, 'CG', problem.Pu)
        Q = FunctionSpace(mesh, 'CG', problem.Pp)
        R = FunctionSpace(mesh, 'CG', problem.Pp)
        W = MixedFunctionSpace([V, Q, R])

        #define trial and test function
        v, q, r = TestFunctions(W)

        w = Function(W)
        w_ = Function(W)

        #initial condition
        w = InitialConditions(problem, V, Q, R) 
        w = project(w, W)


        #initial condition
        w = InitialConditions(problem, V, Q, R) 
        w = project(w, W)

        w_ = InitialConditions(problem, V, Q, R) 
        w_ = project(w_, W)

        U, p, rho = (as_vector((w[0], w[1])), w[2], w[3])
        U_, p_, rho_ = (as_vector((w_[0], w_[1])), w_[2], w_[3])

        # Get boundary conditions
        bcs = problem.boundary_conditions(W.sub(0), W.sub(1), t)

        #rho_(k+theta)
        rho_theta = (1.0-theta)*rho_ + theta*rho

        #U_(k+theta)
        U_theta = (1.0-theta)*U_ + theta*U

        #p_(k+theta)
        p_theta = (1.0-theta)*p_ + theta*p

        f = problem.F
        #weak form of the equations
        #density equation
        F = ((1./dt)*(rho - rho_) + inner(U_theta,grad(rho_theta)))*r*dx
        #momentum equation
        F += (rho_theta/dt)*inner(U - U_,v)*dx \
            + rho_theta*inner(grad(U_theta)*U_theta,v)*dx \
            + nu*inner(grad(U_theta),grad(v))*dx \
            + inner(grad(p_theta),v)*dx 
        #load vector
        F -= rho_theta*inner(theta*f(t) + (1. - theta)*f(t+theta),v)*dx
        #continuity
        F += div(U_theta)*q*dx 
        #stabilization
        if(problem.stabilize):
          # Stabilization parameters
          k1  = 0.5
          k2  = 1.0
          d1 = k1*(dt**(-2) + inner(U_,U_)*h**(-2))**(-0.5) 
          d2 = k2*h 
          #add stabilization
          F += d1*inner(grad(U_theta)*U_theta + grad(p_theta), \
              grad(v)*U_theta + grad(q))*dx + d2*div(U_theta)*div(v)*dx

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
            p_ = w_.split()[1]
            rho_ = w_.split()[2]

            # Update
            self.update(problem, t, U_, rho_)
        return U_, rho_

    def __str__(self):
          return 'DensityNSE'
