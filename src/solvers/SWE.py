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

    #strong residual for cG(1)cG(1)
    def strong_residual(self,u,U,eta):
        H = self.H
        f0 = self.f0
        beta = self.beta
        g = self.g
        nu = self.nu
        rho = self.rho

        NonLinear = self.NonLinear

        R1 = H*div(U)
        R2 = NonLinear*grad(u)*U \
            + f(f0,beta)*as_vector((-U[1],U[0])) \
            + g*grad(eta)

        return R1, R2

    #weak residual for cG(1)cG(1)
    def weak_residual(self,U,U_,eta,eta_,v,chi):
        H = self.H
        f0 = self.f0
        beta = self.beta
        g = self.g
        nu = self.nu
        rho = self.rho

        NonLinear = self.NonLinear

        theta = self.options['theta'] #time stepping method
        dt = self.options["dt"]

        #U_(k+theta)
        U_theta = (1.0-theta)*U_ + theta*U

        #p_(k+theta)
        eta_theta = (1.0-theta)*eta_ + theta*eta

        #weak form of the equations
        r = (1./dt)*(eta - eta_)*chi*dx \
            + H*div(U_theta)*chi*dx 
        r += (1./dt)*inner(U - U_,v)*dx \
            + f(f0,beta)*(U_theta[0]*v[1] - U_theta[1]*v[0])*dx \
            - g*eta_theta*div(v)*dx
        #add the terms for the non-linear SWE
        r += NonLinear*(inner(grad(U_theta)*U_theta,v)*dx \
            + nu*inner(grad(U_theta),grad(v))*dx)

        return r

    def solve(self, problem):
        #get problem mesh 
        mesh = problem.mesh
        h = CellSize(mesh) #mesh size

        #problem parameters
        self.H = problem.h
        self.f0 = problem.f0
        self.beta = problem.beta
        self.g = problem.g
        self.nu = problem.nu
        self.rho = problem.rho

        #if we want a linear version then make a coefficient zero for the
        #terms which only occur in the non-linear from of SWE
        if(self.options['linear']):
            self.NonLinear = 0
        else:
            self.NonLinear = 1

        t = 0 #initial time
        T = problem.T #final time
        dt = self.options['dt'] #time step
        theta = self.options['theta'] #time stepping method
        Pu = self.options["velocity_order"] #order of velocity element
        Pp = self.options["height_order"] #order of height/pressure element

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
        #w = self.InitialConditions(problem, W)
        w_ = self.InitialConditions(problem, W)

        U, eta = (as_vector((w[0], w[1])), w[2])
        U_, eta_ = (as_vector((w_[0], w_[1])), w_[2])

        #U_(k+theta)
        U_theta = (1.0-theta)*U_ + theta*U

        #p_(k+theta)
        eta_theta = (1.0-theta)*eta_ + theta*eta

        F = self.weak_residual(U, U_, eta, eta_, v, chi)

        if(self.options['stabilize']):
          # Stabilization parameters
          k1  = 0.5
          k2  = 0.5
          d1 = k2*(dt**(-2) + eta_*eta_*h**(-2))**(-0.5) 
          d2 = k1*(dt**(-2) + inner(U_,U_)*h**(-2))**(-0.5)

          #add stabilization
          R1, R2 = self.strong_residual(U_theta,U_theta,eta_theta)
          Rv1, Rv2 = self.strong_residual(U_theta,v,chi)
          F += d1*R1*Rv1*dx + d2*inner(R2,Rv2)*dx

        U_, p_ = self.timeStepper(problem, t, T, dt, W, w, w_, U_, eta_, F) 
        return U_, eta_

    def __str__(self):
          return 'SWE'
