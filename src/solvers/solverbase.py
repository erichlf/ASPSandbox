__author__ = "Anders Logg <logg@simula.no>"
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

StabileSolvers = ['NSE', 'SWE', 'MovingSWE']
LinearSolvers = ['SWE']

class SolverBase:
#   Base class for all solvers.
    def __init__(self, options):

        # Store options
        self.options = options

        # Reset some solver variables
        self._time = None
        self._cputime = 0.0
        self._timestep = 0

        # Reset files for storing solution
        self._ufile = None
        self._pfile = None

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
        #get problem mesh 
        mesh = problem.mesh
        h = CellSize(mesh) #mesh size

        #problem parameters
        self.Fr = self.options['Fr'] #Froude Number
        self.Th = self.options['Th'] #Ratio of wave height to depth
        self.Re = self.options['Re'] #Reynolds Number
        self.Ro = self.options['Ro'] #Rossby Number

        #if we want a linear version then make a coefficient zero for the
        #terms which only occur in the non-linear from of SWE
        if(self.options['linear']):
            self.NonLinear = 0
        else:
            self.NonLinear = 1

        t = 0 #initial time
        T = problem.T #final time
        self.dt = self.options['dt'] #time step
        theta = self.options['theta'] #time stepping method
        Pu = self.options['velocity_order'] #order of velocity element
        Pp = self.options['height_order'] #order of height/pressure element

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
          if(self.Fr is None):
              k1  = 1.0
              k2  = 0.5
          else:
              k1  = 1.0*self.Th
              k2  = 0.5*self.Fr**(2.)
          d1 = k1*(dt**(-2) + eta_*eta_*h**(-2))**(-0.5) 
          d2 = k2*(dt**(-2) + inner(U_,U_)*h**(-2))**(-0.5)

          #add stabilization
          R1, R2 = self.strong_residual(U_theta,U_theta,eta_theta)
          Rv1, Rv2 = self.strong_residual(U_theta,v,chi)
          F += d1*R1*Rv1*dx + d2*inner(R2,Rv2)*dx

        U_, eta_ = self.timeStepper(problem, t, T, dt, W, w, w_, U_, eta_, F) 
        return U_, eta_

    def timeStepper(self, problem, t, T, dt, W, w, w_, U_, eta_, F):
        # Time loop
        self.start_timing()

        #plot and save initial condition
        self.update(problem, t, w_.split()[0], w_.split()[1])

        while t<T:
            t += dt

            #evaluate bcs again (in case they are time-dependent)
            bcs = problem.boundary_conditions(W.sub(0), W.sub(1), t)

            solve(F==0, w, bcs=bcs)

            w_.vector()[:] = w.vector()

            U_, eta_ = w_.split()

            # Update
            self.update(problem, t, U_, eta_)

        return U_, eta_

    def prefix(self, problem):
        #Return file prefix for output files
        p = problem.__module__.split('.')[-1]
        s = self.__module__.split('.')[-1]
        if(self.options['linear'] and s in LinearSolvers):
            if(self.options['stabilize'] and s in StabileSolvers):
                s += 'Stabilized'
            s = 'Linear' + s
        elif(self.options['stabilize'] and s in StabileSolvers):
                s += 'Stabilized'

        return problem.output_location + p + s

    def suffix(self):
        #Return file suffix for output files
        if(not self.options['linear']):
            s = 'Re' + str(int(self.Re))
        if(self.Fr is not None):
            s += 'Fr' + str(int(1./self.Fr))
        if(self.Ro is not None):
            s += 'Ro' + str(int(1./self.Ro))
        if(self.Th is not None):
            s += 'Th' + str(int(1./self.Th))

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
            N = self.options['N']
            if (self._timestep - 1) % frequency == 0:
                # Create files for saving
                if self._ufile is None:
                    s = 'results/' + self.prefix(problem) \
                            + self.suffix() \
                            + 'N' + str(N) \
                            + 'K' + str(int(1./dt))
                    self._ufile = File(s + '_u.pvd')
                if self._pfile is None:
                    self._pfile = File(s + '_p.pvd')
                self._ufile << u
                self._pfile << p
        else:
            self.options['plot_solution'] = True

        # Plot solution
        if self.options['plot_solution']:
            if self.vizU is None:
                regex = re.compile('SWE')
                # Plot velocity and pressure
                self.vizU = plot(u, title='Velocity', rescale=True)
                if regex.search(self.prefix(problem)) is None:
                    self.vizP = plot(p, title='Pressure', rescale=True)
                else :
                    self.vizP = plot(p, title='Height', rescale=True)
            else :
                self.vizU.plot(u)
                self.vizP.plot(p,elevate=0.0)

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
        U0, p0 = problem.initial_conditions(W.sub(0),W.sub(1))
        W0 = self.W_project(U0,p0,W)

        return W0
