__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"

from solverbase import *

def f(f0,beta): #Coriolis parameter
    return Expression('f0 + beta*x[1]', f0=f0, beta=beta)

class Solver(SolverBase):
#    Incremental pressure-correction scheme.

    def __init__(self, options):
        SolverBase.__init__(self, options)

    def solve(self, problem):
        #get problem mesh 
        mesh = problem.mesh
        h = CellSize(mesh) #mesh size

        #if we want a linear version then make a coefficient zero for the
        #terms which only occur in the non-linear from of SWE
        if(self.options['linear']):
            linear = 0
        else:
            linear = 1

        t = 0 #initial time
        T = problem.T #final time
        dt = self.options['dt'] #time step
        theta = self.options['theta'] #time stepping method
        Pu = self.options["velocity_order"] #order of velocity element
        Pp = self.options["height_order"] #order of height/pressure element

        #problem parameters
        H = problem.h #fluid depth
        g = problem.g #gravity
        f0 = problem.f0 #reference Coriolis force
        beta = problem.beta #beta plane parameter
        nu = problem.nu #viscosity
        rho = problem.rho #density

        # Define function spaces
        V = VectorFunctionSpace(mesh, 'CG', Pu)
        Q = FunctionSpace(mesh, 'CG', Pp)
        W = MixedFunctionSpace([V, Q])

        # Get boundary conditions
        bcs = problem.boundary_conditions(W.sub(0), W.sub(1), t)

        #define trial and test function
        v, chi = TestFunctions(W)

        w = Function(W)
        w_ = Function(W)

        #initial condition
        w = self.InitialConditions(problem, W)
        w_ = self.InitialConditions(problem, W)

        U, eta = (as_vector((w[0], w[1])), w[2])
        U_, eta_ = (as_vector((w_[0], w_[1])), w_[2])

        #U_(k+theta)
        U_theta = (1.0-theta)*U_ + theta*U

        #p_(k+theta)
        eta_theta = (1.0-theta)*eta_ + theta*eta

        #weak form of the equations
        F = (1./dt)*(eta - eta_)*chi*dx \
            + H*div(U_theta)*chi*dx 
        F += (1./dt)*inner(U - U_,v)*dx \
            + f(f0,beta)*(U_theta[0]*v[1] - U_theta[1]*v[0])*dx \
            - g*eta_theta*div(v)*dx
        #add the terms for the non-linear SWE
        F += linear*(inner(grad(U_theta)*U_theta,v)*dx \
            + nu*inner(grad(U_theta),grad(v))*dx)

        if(self.options['stabilize']):
          # Stabilization parameters
          k1  = 0.5
          k2  = 0.5
          d1 = k1*(dt**(-2) + inner(U_,U_)*h**(-2))**(-0.5)
          d2 = k2*(dt**(-2) + eta_*eta_*h**(-2))**(-0.5) 

          #add stabilization
          F += d1*inner(f(f0,beta)*as_vector((-U_theta[1],U_theta[0])) \
              + g*grad(eta_theta) \
              + linear*grad(U_theta)*U_theta, 
              f(f0,beta)*as_vector((-v[1],v[0])) \
              + g*grad(chi) \
              + linear*grad(v)*U_theta)*dx 
          F += d2*H**2*div(U_theta)*div(v)*dx

        U_, p_ = self.timeStepper(problem, t, T, dt, W, w, w_, U_, eta_, F) 
        return U_, eta_

    def __str__(self):
          return 'SWE'
