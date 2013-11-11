__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"

from solverbase import *

class Solver(SolverBase):
#    Incremental pressure-correction scheme.

    def __init__(self, options):
        SolverBase.__init__(self, options)

    def wave_object(self,W,t):

        """
            Here we create a moving object under water to create a wave.
            Since we are using a mixed space we have to project by solving
            a system. Then we can extract zeta from the solution to the
            system.
        """

        zeta = Expression('A/sqrt(2*S*pi)*exp(-pow(x[0] - (t - 0.5),2)/(2*S))',
                A=200., S=0.02, t=t)

        zeta = self.Q_project(zeta,W)

        return zeta

    def solve(self, problem):
        #get problem mesh
        mesh = problem.mesh
        h = CellSize(mesh) #mesh size

        t = 0 #initial time
        T = problem.T #final time
        dt = problem.dt #time step
        theta = problem.theta #time stepping method

        #parameters
        H = problem.h
        eps = problem.nu
        sigma = problem.nu

        # Define function spaces
        V = VectorFunctionSpace(mesh, 'CG', problem.Pu)
        Q = FunctionSpace(mesh, 'CG', problem.Pp)
        W = MixedFunctionSpace([V, Q])

        #problem parameters
        #H, eps, sigma = get_params(problem, W)

        # Get boundary conditions
        bcs = problem.boundary_conditions(W.sub(0), W.sub(1), t)

        #define trial and test function
        v, chi = TestFunctions(W)

        w = Function(W)
        w_ = Function(W)

        zeta = Function(Q)
        zeta_ = Function(Q)
        zeta__ = Function(Q)

        #initial condition
        w = self.InitialConditions(problem, W))

        w_ = self.InitialConditions(problem, W)

        zeta = self.wave_object(W,t)
        zeta_.assign(zeta)
        zeta__.assign(zeta)

        U, eta = (as_vector((w[0], w[1])), w[2])
        U_, eta_ = (as_vector((w_[0], w_[1])), w_[2])

        #U_(k+theta)
        U_theta = (1.0-theta)*U_ + theta*U

        #p_(k+theta)
        eta_theta = (1.0-theta)*eta_ + theta*eta

        #weak form of the equations
        F = (1./dt)*(eta - eta_)*chi*dx \
            + (H + eps*eta_theta)*div(U_theta)*chi*dx \
            + inner(grad(H + eps*eta_theta),U_theta)*chi*dx
        F += 1./dt*(zeta - zeta_)*chi*dx
        F += (1./dt)*inner(U - U_,v)*dx \
            - eta_theta*div(v)*dx \
            + eps*inner(grad(U_theta)*U_theta,v)*dx \
            + 1./dt*sigma**2*H**2/3.*inner(grad(U - U_),grad(v))*dx
        F -= H/2.*inner(grad(1./dt**2*(zeta - 2*zeta_ + zeta__)),v)*dx

        if(problem.stabilize):
          # Stabilization parameters
          k1  = 0.5
          k2  = 0.5
          d1 = k1*(dt**(-2) + inner(U_,U_)*h**(-2))**(-0.5)
          d2 = k2*(dt**(-2) + eta_*eta_*h**(-2))**(-0.5)

          #add stabilization
          F += d1*(inner(grad(H+eps*eta_theta),U_theta) \
                  + (H + eps*eta_theta)*div(U_theta) + 1./dt*(zeta - zeta_)) \
              *(inner(grad(H+eps*eta_theta),v) + (H + eps*eta_theta)*div(v))*dx
          F += d2*inner(grad(eta_theta) + eps*grad(U_theta)*U_theta \
                  + H/2.*grad(1./dt**2*(zeta - 2*zeta_ + zeta__)), \
              grad(chi) + eps*grad(v)*U_theta)*dx

        # Time loop
        self.start_timing()

        #plot and save initial condition
        self.update(problem, t, w_.split()[0], w_.split()[1])

        while t<T:
            t += dt
            #update the wave generator
            zeta__.assign(zeta_)
            zeta_.assign(zeta)
            zeta = self.wave_object(W,t)

            #evaluate bcs again (in case they are time-dependent)
            bcs = problem.boundary_conditions(W.sub(0), W.sub(1), t)

            solve(F==0, w, bcs=bcs)

            w_.vector()[:] = w.vector()

            U_ = w_.split()[0]
            eta_ = w_.split()[1]

            # Update
            self.update(problem, t, U_, eta_)

        return U_, eta_

    def __str__(self):
          return 'MovingSWE'
