__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"
#
#   adapted from solverbase.py in nsbench originally developed by
#   Anders Logg <logg@simula.no>
#

from dolfin import *

from time import time
from os import getpid
from commands import getoutput
import re
import sys

# Common solver parameters
maxiter = default_maxiter = 200
tolerance = default_tolerance = 1e-4

class SolverBase:
#   Base class for all solvers.
    def __init__(self, options):

        # Store options
        self.options = options

        #initialize parameters
        self.Re = self.options['Re'] #Reynolds number
        self.H = None #Fluid depth
        self.Ro = None #Rossby number
        self.Fr = None #Froude number
        self.Th = None #average wave height
        self.zeta = None #average wave height

        #initialize the time stepping method parameters
        self.t0 = 0 #initial time
        self.dt = self.options['dt'] #time step
        self.alpha = self.options['alpha'] #time stepping method
        self.T = self.options['T'] #Final Time

        #initialize element orders
        self.Pu = self.options['velocity_order'] #order of velocity element
        self.Pp = self.options['height_order'] #order of height/pressure element

        # Reset some solver variables
        self._time = None
        self._cputime = 0.0
        self._timestep = 0

        # Reset files for storing solution
        self._ufile = None
        self._pfile = None
        self._bfile = None

        # Reset storage for functional values and errors
        self._t = []

        #create plot objects
        self.vizU = None
        self.vizP = None

    def getMyMemoryUsage(self):
        mypid = getpid()
        mymemory = getoutput('ps -o rss %s' % mypid).split()[1]
        return mymemory

    def start_timing(self):
#       Start timing, will be paused automatically during update
#       and stopped when the end-time is reached.
        self._time = time()

    def solve(self, problem):
        self.problem = problem
        #get problem mesh
        mesh = problem.mesh
        self.mesh = mesh
        h = CellSize(mesh) #mesh size

        #if we want a linear version then make a coefficient zero for the
        #terms which only occur in the non-linear from of SWE
        if(self.options['linear']):
            self.NonLinear = 0
        else:
            self.NonLinear = 1

        if(self.options['inviscid']):
            self.inviscid = 0
        else:
            self.inviscid = 1

        t = self.t0
        T=self.T
        T = problem.T #final time

        F1 = problem.F1(t) #forcing function for the momentum equation
        F2 = problem.F2(t) #mass source for the continuity equation

        # Define function spaces
        self.V = VectorFunctionSpace(mesh, 'CG', self.Pu)
        V = self.V
        self.Q = FunctionSpace(mesh, 'CG', self.Pp)
        Q = self.Q
        W = MixedFunctionSpace([V, Q])

        # Get boundary conditions
        bcs = problem.boundary_conditions(W.sub(0), W.sub(1), t)

        #forcing and mass source/sink
        F1 = self.V_project(problem.F1(t),W)
        F2 = self.Q_project(problem.F2(t),W)

        #define trial and test function
        v, chi = TestFunctions(W)

        w = Function(W)
        w_ = Function(W)

        #initial condition
        #w_ = self.InitialConditions(problem, W)

        U, eta = (as_vector((w[0], w[1])), w[2])
        U_, eta_ = self.InitialConditions(problem, W)

        #U_(k+alpha)
        U_alpha = (1.0-self.alpha)*U_ + self.alpha*U

        #p_(k+alpha)
        eta_alpha = (1.0-self.alpha)*eta_ + self.alpha*eta

        F = self.weak_residual(U, U_, eta, eta_, v, chi) \
            - inner(F1,v)*dx - F2*chi*dx

        if(self.options['stabilize'] and 'stabilization_parameters' in dir(self)):
          # Stabilization parameters
          d1, d2 = self.stabilization_parameters(U_,eta_,h)

          #add stabilization
          if 'wave_object' in dir(self):
              #R1, R2, z1, z2 = self.strong_residual(U_alpha,U_alpha,eta_alpha)
              #Rv1, Rv2, zv1, zv2 = self.strong_residual(U_alpha,v,chi)
              #F += d1*inner(R1 + z1 - F1, Rv1)*dx + d2*(R2 + z2 - F2)*Rv2*dx
              F += h**(3./2.)*(inner(grad(U_alpha),grad(v)) + inner(grad(eta),grad(chi)))*dx
          else:
              R1, R2 = self.strong_residual(U_alpha,U_alpha,eta_alpha)
              Rv1, Rv2 = self.strong_residual(U_alpha,v,chi)
              F += d1*inner(R1 - F1, Rv1)*dx + d2*(R2 - F2)*Rv2*dx

        U_, eta_ = self.timeStepper(problem, t, T, self.dt, W, w, w_, U_, eta_, F)
        return U_, eta_

    def timeStepper(self, problem, t, T, dt, W, w, w_, U_, eta_, F):
        # Time loop
        self.start_timing()

        #plot and save initial condition
        self.update(problem, t, w_.split()[0], w_.split()[1])

        while t<T:
            t += dt

            if('wave_object' in dir(self)):
                self.wave_object(t, dt)

            #evaluate bcs again (in case they are time-dependent)
            bcs = problem.boundary_conditions(W.sub(0), W.sub(1), t)

            solve(F==0, w, bcs=bcs)

            #update forcing and mass source/sink
            F1 = self.V_project(self.problem.F1(t),W)
            F2 = self.Q_project(self.problem.F2(t),W)

            w_.vector()[:] = w.vector()

            U_, eta_ = w_.split()

            # Update
            self.update(problem, t, U_, eta_)

        return U_, eta_

    def prefix(self, problem):
        #Return file prefix for output files
        p = problem.__module__.split('.')[-1]
        s = self.__module__.split('.')[-1]
        if(self.options['stabilize'] and 'stabilization_parameters' in dir(self)):
            s += 'Stabilized'
        if(self.options['inviscid']):
            s = 'Inviscid' + s
        if(self.options['linear']):
            s = 'Linear' + s

        return problem.output_location + p + s

    def suffix(self):
        s = ''

        #Return file suffix for output files
        if(not self.options['inviscid'] and self.Re is not None):
            s = 'Re' + str(int(self.Re))
        if(self.Ro is not None):
            s += 'Ro' + str(self.Ro)
        if(self.Fr is not None):
            s += 'Fr' + str(self.Fr)
        if(self.Th is not None):
            s += 'Th' + str(self.Th)
        if(self.H is not None):
            s += 'H' + str(self.H)

        return s

    def update(self, problem, t, u, p):
        #Update problem at time t

        # Add to accumulated CPU time
        timestep_cputime = time() - self._time
        self._cputime += timestep_cputime

        # Update problem
        problem.update_problem(t, u, p)

        # Store values
        self._t.append(t)

        # Save solution
        if self.options['save_solution']:
            # Save velocity and pressure
            frequency = self.options['save_frequency']
            dt = self.dt
            N = self.options['Ny']
            if (self._timestep - 1) % frequency == 0:
                # Create files for saving
                if self._ufile is None:
                    s = 'results/' + self.prefix(problem) \
                            + self.suffix() \
                            + 'Ny' + str(N) \
                            + 'K' + str(int(1./dt))
                    self._ufile = File(s + '_u.pvd')
                if self._pfile is None:
                    self._pfile = File(s + '_p.pvd')
                if self._bfile is None:
                    self._bfile = File(s + '_b.pvd')
                self._ufile << u
                self._pfile << p
                self._bfile << self.zz_
        else:
            self.options['plot_solution'] = True

        # Plot solution
        if self.options['plot_solution']:
            if self.vizU is None:
                regex = re.compile('NSE')
                # Plot velocity and pressure
                self.vizU = plot(u, title='Velocity', rescale=True)
                if regex.search(self.prefix(problem)) is None:
                    self.vizP = plot(p, title='Height', rescale=True)
                else :
                    self.vizP = plot(p, title='Pressure', rescale=True, elevate=0.0)
                if('wave_object' in dir(self)):
                    self.vizZ = plot(interpolate(self.h,self.Q), title='Wave Object', rescale=True)
            else :
                self.vizU.plot(u)
                self.vizP.plot(p)
                if('wave_object' in dir(self)):
                    self.vizZ.plot(interpolate(self.h,self.Q))

        # Check memory usage
        if self.options['check_mem_usage']:
            print 'Memory usage is:' , self.getMyMemoryUsage()

        # Print progress
        s = 'Time step %d finished in %g seconds, %g%% done (t = %g, T = %g).' \
            % (self._timestep, timestep_cputime, 100.0*(t / problem.T), t, problem.T)
        sys.stdout.flush()
        sys.stdout.write('\033[K')
        sys.stdout.write(s + '\r')

        # Increase time step and record current time
        self._timestep += 1
        self._time = time()

    def W_project(self,f1,f2,W):
        #This function will project an expressions into W
        e1, e2 = TestFunctions(W)
        w = TestFunction(W)

        u = Function(W)

        A = inner(u,w)*dx - inner(f1,e1)*dx - f2*e2*dx

        solve(A == 0, u)

        return u

    def V_project(self,f1,W):
        #This function will project an expression into W.sub(1)

        #filler function
        f2 = Expression('0.0')

        f1 = self.W_project(f1,f2,W)
        f1 = f1.split()[0]

        return f1

    def Q_project(self,f2,W):
        #This function will project an expression into W.sub(1)

        #filler function
        f1 = Expression(('0.0','0.0'))

        f2 = self.W_project(f1,f2,W)
        f2 = f2.split()[1]

        return f2

    def InitialConditions(self,problem,W):
        #project the given initial condition into W
        U0, p0 = problem.initial_conditions(self.V,self.Q)
        W0 = Function(W)
        U_0, p_0 = W.split() 
        U_0 = interpolate(U0, self.V)
        p_0 = interpolate(p0,self.Q)
        
        return U_0, p_0
