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
from numpy import absolute
import re
import sys
import math

dolfin.parameters["adjoint"]["record_all"] = True

# Common solver parameters
maxiter = default_maxiter = 200
tolerance = default_tolerance = 1e-4

class SolverBase:
    '''
    SolverBase provides a general solver class for the various Solvers. Its
    purpose is take a weak_residual and then solve it using the theta-method
    '''
    def __init__(self, options):

        # Store options
        self.options = options

        #initialize parameters
        self.Re = None #Reynolds number
        self.H = None #Fluid depth
        self.Ro = None #Rossby number
        self.Fr = None #Froude number
        self.Th = None

        #initialize the time stepping method parameters
        self.alpha = self.options['alpha'] #time stepping method

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
        self.s = None
        self._ufile = None
        self._pfile = None
        self._uDualfile = None
        self._pDualfile = None
        self.eifile = None
        self.meshfile = None
        self.optfile = None

        # Reset storage for functional values and errors
        self._t = []

        #create plot objects
        self.vizU = None
        self.vizP = None

    def solve(self, problem):
        '''
            This is the general solve class which will determine if adaptivity
            should be used or if a problem is an optimization problem.
        '''
        mesh = problem.mesh
        k = problem.k #time step
        T = problem.T
        t0 = problem.t0

        TOL = 0*1E-10
        COND = 1

        #naming scheme
        self.s = 'results/' + self.prefix(problem) \
                + self.suffix(problem)

        maxadaps = 2*10 #max number of adaptive steps
        adapt_ratio = 0.1 #number of cells to refine
        nth = ('st','nd','rd','th') #numerical descriptors

        if self.options['adaptive']: #solve with adaptivity
            #files specific to adaptivity
            self.eifile = File(self.s + '_ei.pvd') #error indicators

            # Adaptive loop
            i = 0
            m = 0 #initialize
            while(i<=maxadaps and COND>TOL):
                #setup file names
                self.file_naming(n=i, dual=False)
                #save our current mesh
                if not self.options['plot_solution']:
                    self.meshfile << mesh

                if i==0:
                    print 'Solving on initial mesh.'
                elif i < len(nth):
                    print 'Solving on %d%s adapted mesh.' % (i, nth[i-1])
                else:
                    print 'Solving on %d%s adapted mesh.' % (i, nth[-1])

                # Solve primal and dual problems and compute error indicators
                m_ = m #save the previous functional value
                W, w, m, ei = self.adaptive_solve(problem, mesh, k)
                COND = self.condition(ei, m, m_)
                print 'Stopping Criterion=%0.3G' % COND

                if i==0 and self.options['plot_solution']:
                    plot(ei, title="Error Indicators.", elevate=0.0)
                    plot(mesh, title='Initial mesh', size=((600, 300)))
                elif (i==maxadaps or COND<=TOL) and self.options['plot_solution']:
                    plot(ei, title="Error Indicators.", elevate=0.0)
                    plot(mesh, title='Final mesh', size=((600, 300)))
                    interactive()
                elif not self.options['plot_solution']:
                    self.eifile << ei

                # Refine the mesh
                print 'Refining mesh.'
                mesh = self.adaptive_refine(mesh, ei, adapt_ratio)
                if 'time_step' in dir(self):
                    k = self.time_step(problem.Ubar, mesh)

                adj_reset() #reset the dolfin-adjoint

                i += 1

            if i>maxadaps and COND>TOL:
                s = 'Warning reached max adaptive iterations with sum(abs(EI))=%0.3G.  Solution may not be accurate.' % COND
                print s
        print 'Solving the primal problem.'
        self.file_naming(n=-1, dual=False)

        #record so that we can evaluate our functional
        parameters["adjoint"]["stop_annotating"] = False
        W, w, m = self.forward_solve(problem, mesh, k, func=self.options['adaptive'])
        adj_html("forward.html", "forward")
        if m is not None:
            print
            print 'The size of the functional is: %0.3G' % m

        #solve the optimization problem
        if(self.options['optimize'] and  'Optimize' in dir(self)):
            self.Optimize(problem, w)

        return w

    def adaptive_solve(self, problem, mesh, k):
        '''
            Adaptive solve applies the error representation to goal-oriented
            adaptivity. This is all done automatically using the weak_residual.
        '''
        print 'Solving the primal problem.'
        parameters["adjoint"]["stop_annotating"] = False
        W, w, m = self.forward_solve(problem, mesh, k, func=True)
        parameters["adjoint"]["stop_annotating"] = True
        self._timestep = 0 #reset the time step to zero

        T = problem.T
        t0 = problem.t0

        phi = Function(W)

        # Generate error indicators
        Z = FunctionSpace(mesh, "DG", 0)
        ei = Function(Z, name='Error Indicator')
        z = TestFunction(Z)
        LR1 = 0.

        # Generate the dual problem
        J = Functional(self.functional(mesh, w)*dt[0.9*(T-t0):T], name='DualArgument')
        print
        print 'The size of the functional is: %0.3G' % m
        timestep = None
        wtape = []
        phi = []

        print 'Solving the dual problem.'
        for (adj, var) in compute_adjoint(J,forget=False):
            if var.name == 'w' and timestep != var.timestep:
                timestep = var.timestep
                # Compute error indicators ei
                wtape.append(DolfinAdjointVariable(w).tape_value(timestep=timestep))
                phi.append(adj)
                self.update(problem, None, W, adj, dual=True)

        self._timestep = 0 #reset the time step to zero

        print 'Building error indicators.'
        for i in range(0, len(wtape)-1):
            #the tape is backwards so i+1 is the previous time step
            LR1 = k*self.weak_residual(problem, k, W, wtape[i], wtape[i+1], phi[i], ei_mode=True)
            ei.vector()[:] += assemble(LR1,annotate=False).array()

        return W, w, m, ei

    def condition(self, ei, m, m_):
        '''
            Adaptive stopping criterion for non-Galerkin-orthogonal problems.
            Overload this for Galerkin-orthogonal problems.
            ei - error indicators (non-Galerkin-orthogonal problems)
            m - current functional size (Galerkin-orthogonal problems)
            m_ - previous functional size (Galerkin-orthogonal problems)
        '''
        return abs(sum(ei.vector()))

    def forward_solve(self, problem, mesh, k, func=False):
        '''
            Here we take the weak_residual and apply boundary conditions and then
            send it to time_stepper for solving.
        '''
        h = CellSize(mesh) #mesh size

        t = problem.t0
        T = problem.T #final time

        # Define function spaces
        #we do it this way so that it can be overloaded
        W = self.function_space(mesh)

        # Get boundary conditions
        bcs = problem.boundary_conditions(W, t)

        ic = problem.initial_conditions(W)

        #define trial and test function
        wt = TestFunction(W)
        w = Function(W, name='w')
        w_ = Function(ic, name='w_')

        #weak form of the primal problem
        F = self.weak_residual(problem, k, W, w, w_, wt, ei_mode=False)

        w, m = self.timeStepper(problem, t, T, k, W, w, w_, F, func=func)

        return W, w, m

    #define functions spaces
    def function_space(self, mesh):
        '''
            Sets up a general mixed function space. We assume there are only two
            variables and the first variable is vector valued. To use something
            different the user can overload this in their Solver.
        '''
        V = VectorFunctionSpace(mesh, 'CG', self.Pu)
        Q = FunctionSpace(mesh, 'CG', self.Pp)
        W = MixedFunctionSpace([V, Q])

        return W

    def functional(self, mesh, w):
        '''
            This is the functional used for adaptivity.
            We assume the problem is much like NSE. This can be overloaded by
            each individual problem.
        '''
        if mesh.topology().dim() == 2:
            (u, p) = (as_vector((w[0], w[1])), w[2])
        else:
            (u, p) = (as_vector((w[0], w[1], w[2])), w[3])

        M = u[0]*dx # Mean of the x-velocity in the whole domain

        return M

    # Refine the mesh based on error indicators
    def adaptive_refine(self, mesh, ei, adapt_ratio):
        '''
            Take a mesh and the associated error indicators and refine adapt_ration%
            of cells.
        '''
        gamma = abs(ei.vector().array())

        # Mark cells for refinement
        cell_markers = MeshFunction("bool", mesh, mesh.topology().dim())
        gamma_0 = sorted(gamma, reverse=True)[int(len(gamma)*adapt_ratio) - 1]
        for c in cells(mesh):
            cell_markers[c] = gamma[c.index()] > gamma_0

        # Refine mesh
        mesh = refine(mesh, cell_markers)

        return mesh

    def timeStepper(self, problem, t, T, k, W, w, w_, F, func=False):
        '''
            Time stepper for solver using theta-method.
        '''
        if func:
            m = 0
        else:
            m = None
        bcs = problem.boundary_conditions(W, t)
        # Time loop
        self.start_timing()

        #plot and save initial condition
        self.update(problem, t, W, w_)

        adj_start_timestep(t)
        while t<(T-k/2.):
            t += k

            if('wave_object' in dir(self)):
                self.wave_object(problem, self.Q, t, k)

            #evaluate bcs again (in case they are time-dependent)
            bcs = problem.boundary_conditions(W, t)

            solve(F==0, w, bcs=bcs)

            w_.assign(w)
            if func and t>0.9*T:
                m += 1./k*assemble(self.functional(W.mesh(), w_), annotate=False)
            adj_inc_timestep(t,finished=t>=(T-k/2.))

            # Update
            self.update(problem, t, W, w_)

        return w_, m

    def update(self, problem, t, W, w, dual=False):
        '''
            Saves or plots the data at each time step.
        '''
        # Add to accumulated CPU time
        timestep_cputime = time() - self._time
        self._cputime += timestep_cputime

        # Store values
        if t is not None:
            self._t.append(t)

        # Save solution
        if self.options['save_solution']:
            # Save velocity and pressure
            self.Save(problem, w, dual=dual)
        else: # Plot solution
            self.options['plot_solution'] = True
            self.Plot(problem, W, w)

        # Check memory usage
        if self.options['check_mem_usage']:
            print 'Memory usage is:' , self.getMyMemoryUsage()

        # Print progress
        if not dual:
            s = 'Time step %d finished in %g seconds, %g%% done (t = %g, T = %g).' \
                % (self._timestep, timestep_cputime, 100.0*(t / problem.T), t, problem.T)
            sys.stdout.flush()
            sys.stdout.write('\033[K')
            sys.stdout.write(s + '\r')

            #record current time
            self._time = time()

        # Increase time step
        self._timestep += 1

    def Save(self, problem, w, dual=False):
        '''
            Save a variables associated with a time step. Here we assume there are
            two variables where the first variable is vector-valued and the second
            variable is a scalar. If this doesn't fit the particular solvers
            variables the user will need to overload this function.
        '''
        u = w.split()[0]
        p = w.split()[1]

        if self.options['save_frequency'] != 0 and (self._timestep - 1) % self.options['save_frequency'] == 0:
            if not dual:
                self._ufile << u
                self._pfile << p
            else:
                self._uDualfile << u
                self._pDualfile << p

    def file_naming(self, n=-1, dual=False):
        '''
            Names our files for saving variables. Here we assume there are
            two variables where the first variable is vector-valued and the second
            variable is a scalar. If this doesn't fit the particular solvers
            variables the user will need to overload this function.
        '''
        if n==-1:
            self._ufile = File(self.s + '_u.pvd', 'compressed')
            self._pfile = File(self.s + '_p.pvd', 'compressed')
            self._uDualfile = File(self.s + '_uDual.pvd', 'compressed')
            self._pDualfile = File(self.s + '_pDual.pvd', 'compressed')
            self.meshfile = File(self.s + '_mesh.xml')
        else:
            self._ufile = File(self.s + '_u%d.pvd' % n, 'compressed')
            self._pfile = File(self.s + '_p%d.pvd' % n, 'compressed')
            self._uDualfile = File(self.s + '_uDual%d.pvd' % n, 'compressed')
            self._pDualfile = File(self.s + '_pDual%d.pvd' % n, 'compressed')
            self.meshfile = File(self.s + '_mesh%d.xml' % n)

    #this is a separate function so that it can be overloaded
    def Plot(self, problem, W, w):
        '''
            Plots our variables associated with a time step. Here we assume there are
            two variables where the first variable is vector-valued and the second
            variable is a scalar. If this doesn't fit the particular solvers
            variables the user will need to overload this function.
        '''
        u = w.split()[0]
        p = w.split()[1]

        if self.vizU is None:
            # Plot velocity and pressure
            self.vizU = plot(u, title='Velocity', rescale=True)
            self.vizP = plot(p, title='Pressure', rescale=True, elevate=0.0)
        else :
            self.vizU.plot(u)
            self.vizP.plot(p)

    def prefix(self, problem):
        '''
            Obtains the beginning of file naming, e.g. Probem Name, Solver Name,
            dimension, etc.
        '''
        #Return file prefix for output files
        p = problem.__module__.split('.')[-1]
        if problem.mesh.topology().dim() > 2:
            p += '3D'
        else:
            p += '2D'
        s = self.__module__.split('.')[-1]
        if(self.options['stabilize'] and 'stabilization_parameters' in dir(self)):
            s += 'Stabilized'
        if(self.options['inviscid']):
            s = 'Inviscid' + s
        if(self.options['linear']):
            s = 'Linear' + s

        return problem.output_location + s + p

    def suffix(self, problem):
        '''
            Obtains the run specific data for file naming, e.g. Nx, k, etc.
        '''
        s = ''

        #Return file suffix for output files
        if not self.options['inviscid'] and self.Re is not None:
            s = 'Re' + str(int(self.Re))
        if self.Ro is not None:
            s += 'Ro' + str(self.Ro)
        if self.Fr is not None:
            s += 'Fr' + str(self.Fr)
        if self.Th is not None:
            s += 'Th' + str(self.Th)

        s += 'Nx' + str(problem.Nx)
        if problem.Ny is not None:
            s += 'Ny' + str(problem.Ny)
        if problem.mesh.topology().dim() > 2 and problem.Nz is not None:
            s += 'Nz' + str(problem.Nz)

        s += 'K' + str(int(1./problem.k))

        return s

    def getMyMemoryUsage(self):
        '''
            Determines how much memory we are using.
        '''
        mypid = getpid()
        mymemory = getoutput('ps -o rss %s' % mypid).split()[1]
        return mymemory

    def start_timing(self):
        '''
            Start timing, will be paused automatically during update
            and stopped when the end-time is reached.
        '''
        self._time = time()
