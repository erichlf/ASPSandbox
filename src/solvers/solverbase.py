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
        mesh = problem.mesh
        k = problem.k #time step
        T = problem.T
        t0 = problem.t0

        #naming scheme
        self.s = 'results/' + self.prefix(problem) \
                + self.suffix(problem)

        maxadaps = 10 #max number of adaptive steps
        adapt_ratio = 0.1 #number of cells to refine
        nth = ('st','nd','rd','th') #numerical descriptors

        if self.options['adaptive']: #solve with adaptivity
            #files specific to adaptivity
            self.eifile = File(self.s + '_ei.pvd') #error indicators
            self.meshfile = File(self.s + '_mesh.pvd')

            # Adaptive loop
            for i in range(0, maxadaps):
                #setup file names
                self.file_naming(n=i, dual=False)
                if i==0:
                    print 'Solving on initial mesh.'
                elif i < len(nth):
                    print 'Solving on %d%s adapted mesh.' % (i, nth[i-1])
                else:
                    print 'Solving on %d%s adapted mesh.' % (i, nth[-1])
                # Solve primal and dual problems and compute error indicators
                w, ei = self.adaptive_solve(problem, mesh, k)
                if i==0 and self.options['plot_solution']:
                    plot(ei, title="Error Indicators.", elevate=0.0)
                    plot(mesh, title='Initial mesh', size=((600, 300)))
                elif i==maxadaps-1 and self.options['plot_solution']:
                    plot(ei, title="Error Indicators.", elevate=0.0)
                    plot(mesh, title='Final mesh', size=((600, 300)))
                    interactive()
                elif not self.options['plot_solution']:
                    self.meshfile << mesh
                    self.eifile << ei

                # Refine the mesh
                print 'Refining mesh.'
                mesh = self.adaptive_refine(mesh, ei, adapt_ratio)
                self._timestep = 0 #reset the time step to zero
                adj_reset() #reset the dolfin-adjoint

        print 'Solving the primal problem.'
        self.file_naming(n=-1, dual=False)

        #record so that we can evaluate our functional
        parameters["adjoint"]["stop_annotating"] = False
        W, w = self.forward_solve(problem, mesh, k)
        J = Functional(self.functional(mesh, w)*dt[0.9*(T-t0):T], name='DualArgument')
        m = ScalarParameter('w_')
        RJ = ReducedFunctional(J, m)
        print
        print 'The size of the functional is: %f' % RJ(w)

        #solve the optimization problem
        if(self.options['optimize'] and  'Optimize' in dir(self)):
            self.optfile = File(self.s + '_Opt.pvd') #file for solution to optimization

            opt = self.Optimize(problem, w)
            if self.options['plot_solution']:
                plot(opt, title='Optimization result.')
                interactive()
            else:
                optfile << opt

        return w

    def adaptive_solve(self, problem, mesh, k):

        print 'Solving the primal problem.'
        parameters["adjoint"]["stop_annotating"] = False
        W, w = self.forward_solve(problem, mesh, k)
        parameters["adjoint"]["stop_annotating"] = True

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
        m = ScalarParameter('w_')
        RJ = ReducedFunctional(J, m)
        print
        print 'The size of the functional is: %f' % RJ(w)
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
                if self.options['save_solution']:
                    self.update(problem, None, W, adj, dual=True)

        print 'Building error indicators.'
        for i in range(0, len(wtape)-1):
            LR1 = k*self.weak_residual(problem, W, wtape[i], wtape[i+1], phi[i], ei_mode=True)
            ei.vector()[:] += assemble(LR1,annotate=False).array()

        return w, ei

    def forward_solve(self, problem, mesh, k):
        h = CellSize(mesh) #mesh size

        t = problem.t0
        T = problem.T #final time
        self.k = k

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
        F = self.weak_residual(problem, W, w, w_, wt, ei_mode=False)

        w = self.timeStepper(problem, t, T, k, W, w, w_, F)

        return W, w

    #define functions spaces
    def function_space(self, mesh):
        V = VectorFunctionSpace(mesh, 'CG', self.Pu)
        Q = FunctionSpace(mesh, 'CG', self.Pp)
        W = MixedFunctionSpace([V, Q])

        return W

    def InitialConditions(self,problem,W):
        #project the given initial condition into W
        U0, p0 = problem.initial_conditions(W.sub(0),W.sub(1))
        W0 = self.W_project(U0,p0,W)

        return W0

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

            #update forcing and mass source/sink
            #F1 = self.V_project(self.problem.F1(t),W)
            #F2 = self.Q_project(self.problem.F2(t),W)

            w_.assign(w)
            adj_inc_timestep(t,finished=t>=(T-k/2.))

            # Update
            self.update(problem, t, W, w_)

        return w_

    def update(self, problem, t, W, w, dual=False):
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

            # Increase time step and record current time
            self._timestep += 1
            self._time = time()
#        else:
#            s = 'Calculating the dual finished in %g seconds.' % timestep_cputime
#            sys.stdout.flush()
#            sys.stdout.write('\033[K')
#            sys.stdout.write(s + '\r')
#
#            # Increase time step and record current time
#            self._timestep += 1
#            self._time = time()

    def Save(self, problem, w, dual=False):
        u = w.split()[0]
        p = w.split()[1]

        if (self._timestep - 1) % self.options['save_frequency'] == 0:
            if not dual:
                self._ufile << u
                self._pfile << p
            else:
                self._uDualfile << u
                self._pDualfile << p

    def file_naming(self, n=-1, dual=False):
        if n==-1:
            self._ufile = File(self.s + '_u.pvd', 'compressed')
            self._pfile = File(self.s + '_p.pvd', 'compressed')
            self._uDualfile = File(self.s + '_uDual.pvd', 'compressed')
            self._pDualfile = File(self.s + '_pDual.pvd', 'compressed')
        else:
            self._ufile = File(self.s + '_u%d.pvd' % n, 'compressed')
            self._pfile = File(self.s + '_p%d.pvd' % n, 'compressed')
            self._uDualfile = File(self.s + '_uDual%d.pvd' % n, 'compressed')
            self._pDualfile = File(self.s + '_pDual%d.pvd' % n, 'compressed')

    #this is a separate function so that it can be overloaded
    def Plot(self, problem, W, w):
        u = w.split()[0]
        p = w.split()[1]

        if self.vizU is None:
            # Plot velocity and pressure
            self.vizU = plot(u, title='Velocity', rescale=True)
            self.vizP = plot(p, title='Height', rescale=True, elevate=0.0)
        else :
            self.vizU.plot(u)
            self.vizP.plot(p)

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

    def suffix(self, problem):
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

        s += 'Nx' + str(problem.Nx) \
            + 'Ny' + str(problem.Ny)
        if problem.mesh.topology().dim() > 2:
            s += 'Nz' + str(problem.Nz) \

        s += 'K' + str(int(1./problem.k))

        return s

    def getMyMemoryUsage(self):
        mypid = getpid()
        mymemory = getoutput('ps -o rss %s' % mypid).split()[1]
        return mymemory

    def start_timing(self):
#       Start timing, will be paused automatically during update
#       and stopped when the end-time is reached.
        self._time = time()
