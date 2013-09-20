__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"

from solverbase import *

class InitialConditions(Expression):
    def __init__(self, problem, V, Q, A):
        C = 1
        sig = 0.05
        self.U0, self.eta0 = problem.initial_conditions(V, Q)
        self.U0 = project(self.U0,V)
        self.p0 = project(self.eta0,Q)
        #self.a0 = Expression('C*exp(-(pow(x[0]-0.75,2)+pow(x[1]-0.75,2))/(2*sig*sig))', C=C, sig=sig)
        self.a0 = Expression('(x[0]<0.5)?1.0:0.0')
        self.a0 = project(self.a0,A)

    def eval(self, value, x):
        value[:2] = self.U0(x)
        value[2] = self.p0(x)
        value[3] = self.a0(x)

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
        A = FunctionSpace(mesh, 'CG', problem.Pp)
        W = MixedFunctionSpace([V, Q, A])

        #define trial and test function
        v, q, b = TestFunctions(W)

        w = Function(W)
        w_ = Function(W)

        #initial condition
        w = InitialConditions(problem, V, Q, A) 
        w = project(w, W)


        #initial condition
        w = InitialConditions(problem, V, Q, A) 
        w = project(w, W)

        w_ = InitialConditions(problem, V, Q, A) 
        w_ = project(w_, W)

        U, p, a = (as_vector((w[0], w[1])), w[2], w[3])
        U_, p_, a_ = (as_vector((w_[0], w_[1])), w_[2], w_[3])

        # Get boundary conditions
        bcs = problem.boundary_conditions(W.sub(0), W.sub(1), t)

        #a_(k+theta)
        a_theta = (1.0-theta)*a_ + theta*a

        #U_(k+theta)
        U_theta = (1.0-theta)*U_ + theta*U

        #p_(k+theta)
        p_theta = (1.0-theta)*p_ + theta*p

        f = problem.F
        #weak form of the equations
        F = ((1./dt)*(a - a_) + inner(U_theta,grad(a_theta)))*b*dx
        F += (1./dt)*inner(U - U_,v)*dx \
            + inner(grad(U_theta)*U_theta,v)*dx \
            + nu*(inner(grad(U_theta)*grad(a_theta),v) \
            + (1+a_theta)*inner(grad(U_theta),grad(v)))*dx \
            + (1+a_theta)*inner(grad(p_theta),v)*dx 
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
            a_ = w_.split()[2]

            # Update
            self.update(problem, t, U_, a_)
        return U_, a_

    def __str__(self):
          return 'DensityNSE'
