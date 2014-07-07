__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"
#
#   adapted from solverbase.py in nsbench originally developed by
#   Anders Logg <logg@simula.no>
#

from dolfin import *
from dolfin_adjoint import *

from time import time
from os import getpid
from commands import getoutput
import re
import sys
import math

dolfin.parameters["adjoint"]["record_all"] = True

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
        self.Th = None 
        self.zeta = None

        #initialize the time stepping method parameters
        self.t0 = 0 #initial time
        self.k = self.options['dt'] #time step
        self.alpha = self.options['alpha'] #time stepping method
        self.T = self.options['T'] #Final Time

        #initialize element orders
        self.Pu = self.options['velocity_order'] #order of velocity element
        self.Pp = self.options['height_order'] #order of height/pressure element

        if(self.Pu == 2 and self.options['stabilize']):
            self.options['stabilize'] = False

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
        mesh = problem.mesh

        maxiters = 10 #max number of adaptive steps
        adapt_ratio = 0.1 #number of cells to refine
        nth = ('st','nd','rd','th') #numerical descriptors

        if not self.options['adaptive']: #solve without adaptivity
            u, p = self.forward_solve(mesh)
        else:
            # Adaptive loop
            for i in range(0, maxiters):
                print
                if i==0:
                    print 'Solving on initial mesh.'
                else:
                    if i < len(nth):
                        print 'Solving on %d%s adapted mesh.' % (i, nth[i-1])
                    else:
                        print 'Solving on %d%s adapted mesh.' % (i, nth[len(nth)-1])
                # Solve primal and dual problems and compute error indicators
                (U_, eta_, ei) = self.adaptive_solve(mesh)
                if(i == 0 and self.options['plot_solution']):
                    plot(mesh, title="Initial mesh", size=((600, 300)))
                elif(i == maxiters - 1 and self.options['plot_solution']):
                    plot(mesh, title="Finest mesh", size=((600, 300)))

                # Refine the mesh
                mesh = self.adaptive_refine(mesh, ei, adapt_ratio)
                self._timestep = 0 #reset the time step to zero
                adj_reset() #reset the dolfin-adjoint
            interactive()

        return U_, eta_

    def forward_solve(self,mesh):
        problem = self.problem
        h = CellSize(mesh) #mesh size

        t = self.t0
        T = problem.T #final time
        k = self.k #time step

        # Define function spaces
        V = VectorFunctionSpace(mesh, 'CG', self.Pu)
        Q = FunctionSpace(mesh, 'CG', self.Pp)
        W = MixedFunctionSpace([V, Q])

        # Get boundary conditions
        bcs = problem.boundary_conditions(W.sub(0), W.sub(1), t)

        #define trial and test function
        wt = TestFunction(W)
        v, chi = TestFunctions(W)

        w = Function(W)
        w_ = Function(W)

        #initial condition
        w_ = self.InitialConditions(problem, W)

        U, eta = (as_vector((w[0], w[1])), w[2])
        U_, eta_ = (as_vector((w_[0], w_[1])), w_[2])

        #weak form of the primal problem
        F = self.weak_residual(W, w, w_, wt, ei_mode=False)

        w_ = self.timeStepper(problem, t, T, k, W, w, w_, F)

        return U_, eta_

    def adaptive_solve(self, mesh):
        problem = self.problem
        h = CellSize(mesh) #mesh size

        t = self.t0
        T = problem.T #final time
        k = self.k

        # Define function spaces
        Z = FunctionSpace(mesh, "DG", 0)
        V = VectorFunctionSpace(mesh, 'CG', self.Pu)
        Q = FunctionSpace(mesh, 'CG', self.Pp)
        W = MixedFunctionSpace([V, Q])

        # Get boundary conditions
        bcs = problem.boundary_conditions(W.sub(0), W.sub(1), t)

        #define trial and test function
        wt = TestFunction(W)
        v, chi = TestFunctions(W)

        w = Function(W, name='State')
        w_ = Function(W, name='PreviousState')

        #initial condition
        w_ = self.InitialConditions(problem, W)

        U, eta = (as_vector((w[0], w[1])), w[2])
        U_, eta_ = (as_vector((w_[0], w_[1])), w_[2])

        #weak form of the primal problem
        F = self.weak_residual(W, w, w_, wt, ei_mode=False)

        w_ = self.timeStepper(problem, t, T, k, W, w, w_, F)

        phi = Function(W)

        # Generate error indicators
        ei = Function(Z)
        z = TestFunction(Z)
        LR1 = 0.

        # Generate the dual problem
        J = Functional(self.Functional(mesh, v, chi)*dt)
        i = int(math.ceil(T/k)) #last time step
        adjoint = compute_adjoint(J,forget=False)
        for (phi, var) in adjoint:
          if var.name == 'State':
            # Compute error indicators ei
            wtape = DolfinAdjointVariable(w).tape_value(iteration=i)
            if i>0:
                wtape_ = DolfinAdjointVariable(w).tape_value(iteration=i-1)
            LR1 = k*self.weak_residual(W, wtape, wtape_, phi, ei_mode=True)
            ei.vector()[:] += assemble(LR1).array()
            i -= 1

        return U_, eta_, ei

    # Refine the mesh based on error indicators
    def adaptive_refine(self, mesh, ei, adapt_ratio):
        gamma = abs(ei.vector().array())

        # Mark cells for refinement
        cell_markers = MeshFunction("bool", mesh, mesh.topology().dim())
        gamma_0 = sorted(gamma, reverse=True)[int(len(gamma)*adapt_ratio) - 1]
        for c in cells(mesh):
            cell_markers[c] = gamma[c.index()] > gamma_0

        # Refine mesh
        mesh = refine(mesh, cell_markers)

        return mesh

    def timeStepper(self, problem, t, T, k, W, w, w_, F):
        # Time loop
        self.start_timing()

        #plot and save initial condition
        self.update(problem, t, w_.split()[0], w_.split()[1])

        while t<(T-k/2.):
            t += k

            if('wave_object' in dir(self)):
                self.wave_object(W.sub(1), t, k)

            #evaluate bcs again (in case they are time-dependent)
            bcs = problem.boundary_conditions(W.sub(0), W.sub(1), t)

            solve(F==0, w, bcs=bcs)

            #update forcing and mass source/sink
            F1 = self.V_project(self.problem.F1(t),W)
            F2 = self.Q_project(self.problem.F2(t),W)

            w_.assign(w)

            # Update
            self.update(problem, t, w.split()[0], w.split()[1])

        return w_

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
        #problem.update_problem(t, u, p)

        # Store values
        self._t.append(t)

        # Save solution
        if self.options['save_solution']:
            # Save velocity and pressure
            frequency = self.options['save_frequency']
            k = self.k
            Nx = self.options['Nx']
            Ny = self.options['Ny']
            if (self._timestep - 1) % frequency == 0:
                # Create files for saving
                if self._ufile is None:
                    s = 'results/' + self.prefix(problem) \
                            + self.suffix() \
                            + 'Nx' + str(Nx) \
                            + 'Ny' + str(Ny) \
                            + 'K' + str(int(1./k))
                    self._ufile = File(s + '_u.pvd')
                if self._pfile is None:
                    self._pfile = File(s + '_p.pvd')
                if self._bfile is None:
                    self._bfile = File(s + '_b.pvd')
                self._ufile << u
                self._pfile << p
                self._bfile << self.H_
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
                    self.vizZ = plot(self.H, title='Wave Object', rescale=True)
            else :
                self.vizU.plot(u)
                self.vizP.plot(p)
                if('wave_object' in dir(self)):
                    self.vizZ.plot(self.H_)

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
