__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__ = "GNU GPL version 3 or any later version"
#
#   adapted from solverbase.py in nsbench originally developed by
#   Anders Logg <logg@simula.no>
#

from dolfin import *
from dolfin_adjoint import *

from time import time
from os import getpid
from commands import getoutput
import sys

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

        # initialize parameters
        self.Re = None  # Reynolds number
        self.H = None  # Fluid depth
        self.Ro = None  # Rossby number
        self.Fr = None  # Froude number
        self.Th = None

        # initialize the time stepping method parameters
        self.theta = self.options['theta']  # time stepping method

        # initialize element orders
        self.Pu = self.options['velocity_order']  # order of velocity element
        # order of height/pressure element
        self.Pp = self.options['height_order']

        if(self.Pu == 2 and self.options['stabilize']):
            self.options['stabilize'] = False

        # Reset some solver variables
        self._time = None
        self._cputime = 0.0
        self._timestep = 0

        self.adapt_ratio = self.options['adapt_ratio']
        self.maxadapts = self.options['max_adaptations']
        self.adaptTOL = self.options['adaptive_TOL']

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

        # create plot objects
        self.vizU = None
        self.vizP = None

    def solve(self, problem):
        '''
            This is the general solve class which will determine if adaptivity
            should be used or if a problem is an optimization problem.
        '''
        mesh = problem.mesh
        T = problem.T
        t0 = problem.t0
        # adjust time step so that we evenly divide time interval
        k = self.adjust_dt(t0, T, problem.k)

        # naming scheme
        self.s = 'results/' + self.prefix(problem) + self.suffix(problem)

        if self.options['adaptive']:  # solve with adaptivity
            mesh, k = self.adaptivity(problem, mesh,
                                      T, t0, k)

        print 'Solving the primal problem.'
        self.file_naming(n=-1, dual=False)

        # record so that we can evaluate our functional
        parameters["adjoint"]["stop_annotating"] = \
            not (self.options['adaptive'] or (self.options['optimize']
                 and 'Optimize' in dir(self)))
        W, w, m = self.forward_solve(problem, mesh, t0, T, k,
                                     func=True)#self.options['adaptive'])
        if m is not None:
            print
            print 'The size of the functional is: %0.3G' % m

        # solve the optimization problem
        if(self.options['optimize'] and 'Optimize' in dir(self)):
            self.Optimize(problem, w)

        return w

    def adaptivity(self, problem, mesh, T, t0, k):
        COND = 1
        # files specific to adaptivity
        self.eifile = File(self.s + '_ei.pvd')  # error indicators
        nth = ('st', 'nd', 'rd', 'th')  # numerical descriptors

        # Adaptive loop
        i = 0
        m = 0  # initialize
        while(i <= self.maxadapts and COND > self.adaptTOL):
            # setup file names
            self.file_naming(n=i, dual=False)
            # save our current mesh
            if not self.options['plot_solution']:
                self.meshfile << mesh

            if i == 0:
                print 'Solving on initial mesh.'
            elif i < len(nth):
                print 'Solving on %d%s adapted mesh.' % (i, nth[i - 1])
            else:
                print 'Solving on %d%s adapted mesh.' % (i, nth[-1])

            # Solve primal and dual problems and compute error indicators
            m_ = m  # save the previous functional value
            W, w, m, ei = self.adaptive_solve(problem, mesh, t0, T, k)
            COND = self.condition(ei, m, m_)
            print 'DOFs=%d functional=%0.5G err_est=%0.5G' \
                % (mesh.num_vertices(), m, COND)

            if i == 0 and self.options['plot_solution']:
                plot(ei, title="Error Indicators.", elevate=0.0)
                plot(mesh, title='Initial mesh', size=((600, 300)))
            elif (i == self.maxadapts or COND <= self.adaptTOL) \
                    and self.options['plot_solution']:
                plot(ei, title="Error Indicators.", elevate=0.0)
                plot(mesh, title='Final mesh', size=((600, 300)))
                interactive()
            elif not self.options['plot_solution']:
                self.eifile << ei

            # Refine the mesh
            print 'Refining mesh.'
            mesh = self.adaptive_refine(mesh, ei)
            if 'time_step' in dir(problem):
                k = self.adjust_dt(t0, T, problem.time_step(problem.Ubar, mesh))

            adj_reset()  # reset the dolfin-adjoint

            i += 1

        if i > self.maxadapts and COND > self.adaptTOL:
            s = 'Warning reached max adaptive iterations with' \
                + 'sum(abs(EI))=%0.3G.  Solution may not be accurate.' \
                % COND
            print s

        return mesh, k

    def adaptive_solve(self, problem, mesh, t0, T, k):
        '''
            Adaptive solve applies the error representation to goal-oriented
            adaptivity. This is all done automatically using the weak_residual.
        '''
        print 'Solving the primal problem.'
        parameters["adjoint"]["stop_annotating"] = False

        N = int(round((T - t0) / k))
        perN = self.options['onDisk']

        assert perN <= 1. or perN >= 0.
        if perN > 0:
            adj_checkpointing(strategy='multistage', steps=N,
                              snaps_on_disk=int(perN * N),
                              snaps_in_ram=int((1. - perN) * N),
                              verbose=False)

        W, w, m = self.forward_solve(problem, mesh, t0, T, k, func=True)
        parameters["adjoint"]["stop_annotating"] = True
        self._timestep = 0  # reset the time step to zero

        phi = Function(W)

        # Generate error indicators
        Z = FunctionSpace(mesh, "DG", 0)
        ei = Function(Z, name='Error Indicator')
        LR1 = 0.

        # Generate the dual problem
        J = Functional(self.functional(problem, mesh, w) * dt,
                       name='DualArgument')
        timestep = None
        wtape = []
        phi = []

        print
        print 'Solving the dual problem.'
        for (adj, var) in compute_adjoint(J, forget=False):
            if var.name == 'w' and timestep != var.timestep:
                timestep = var.timestep
                # Compute error indicators ei
                wtape.append(DolfinAdjointVariable(w).
                             tape_value(timestep=timestep))
                phi.append(adj)
                self.update(problem, None, W, adj, dual=True)

        self._timestep = 0  # reset the time step to zero

        print 'Building error indicators.'
        for i in range(0, len(wtape) - 1):
            # the tape is backwards so i+1 is the previous time step
            wtape_theta = self.theta * \
                wtape[i] + (1. - self.theta) * wtape[i + 1]
            LR1 = k * self.weak_residual(problem, k, W, wtape_theta, wtape[i],
                                         wtape[i + 1], phi[i], ei_mode=True)
            ei.vector()[:] += assemble(LR1, annotate=False).array()

        return W, w, m, ei

    def condition(self, ei, m, m_):
        '''
            Adaptive stopping criterion for non-Galerkin-orthogonal problems.
            Overload this for Galerkin-orthogonal problems.
            ei - error indicators (non-Galerkin-orthogonal problems)
            m - current functional size (Galerkin-orthogonal problems)
            m_ - previous functional size (Galerkin-orthogonal problems)
        '''
        if(self.options['stabilize']
           and 'stabilization_parameters' in dir(self)):
            c = abs(sum(ei.vector()))
        else:
            c = abs(m - m_)

        return c

    def forward_solve(self, problem, mesh, t0, T, k, func=False):
        '''
            Here we take the weak_residual and apply boundary conditions and
            then send it to time_stepper for solving.
        '''
        # Define function spaces
        # we do it this way so that it can be overloaded
        W = self.function_space(mesh)

        ic = problem.initial_conditions(W)

        # define trial and test function
        wt = TestFunction(W)
        w = Function(W, name='w')
        w_ = Function(ic, name='w_')

        theta = self.theta
        w_theta = (1. - theta) * w_ + theta * w

        # weak form of the primal problem
        F = self.weak_residual(problem, k, W, w_theta, w, w_, wt, ei_mode=False)

        w, m = self.timeStepper(problem, t0, T, k, W, w, w_, F, func=func)

        return W, w, m

    # define functions spaces
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

    def functional(self, problem, mesh, w):
        '''
            This is the functional used for adaptivity.
            We assume the problem is much like NSE. This can be overloaded by
            each individual problem.
        '''
        if mesh.topology().dim() == 2:
            (u, p) = (as_vector((w[0], w[1])), w[2])
        else:
            (u, p) = (as_vector((w[0], w[1], w[2])), w[3])

        # n = FacetNormal(mesh)
        # marker = problem.marker()

        M = u[0] * dx  # Mean of the x-velocity in the whole domain
        # M = marker*p*n[0]*ds  # Drag (only pressure)

        return M

    # Refine the mesh based on error indicators
    def adaptive_refine(self, mesh, ei):
        '''
            Take a mesh and the associated error indicators and refine
            adapt_ratio% of cells.
        '''
        gamma = abs(ei.vector().array())

        # Mark cells for refinement
        cell_markers = MeshFunction("bool", mesh, mesh.topology().dim())
        gamma_0 = sorted(gamma, reverse=True)[int(len(gamma) * self.adapt_ratio) - 1]
        for c in cells(mesh):
            cell_markers[c] = gamma[c.index()] > gamma_0

        # Refine mesh
        mesh = refine(mesh, cell_markers)

        return mesh

    def adjust_dt(self, t0, T, k):
        '''
            Adjust time step so that we evenly divide the time interval, but
            ensure that the new time step is always smaller than the original.
        '''
        div, rem = divmod((T - t0), k)
        if rem is not 0:
            k = (T - t0) / (div + 1)

        return k

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

        # plot and save initial condition
        self.update(problem, t, W, w_)

        adj_start_timestep(t)
        while t < (T - k / 2.):
            t += k

            if('wave_object' in dir(self)):
                self.wave_object(problem, self.Q, t, k)

            # evaluate bcs again (in case they are time-dependent)
            bcs = problem.boundary_conditions(W, t)

            solve(F == 0, w, bcs=bcs)

            w_.assign(w)
            if func:
                m += k * assemble(self.functional(problem, W.mesh(), w_),
                                  annotate=False)
            adj_inc_timestep(t, finished=t >= (T - k / 2.))

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
        else:  # Plot solution
            self.options['plot_solution'] = True
            self.Plot(problem, W, w)

        # Check memory usage
        if self.options['check_mem_usage']:
            print 'Memory usage is:', self.getMyMemoryUsage()

        # Print progress
        if not dual:
            s = 'Time step %d finished in %g seconds, ' \
                % (self._timestep, timestep_cputime)
            s += '%g%% done (t = %g, T = %g).' % (100.0 * (t / problem.T), t,
                                                  problem.T)
            sys.stdout.flush()
            sys.stdout.write('\033[K')
            sys.stdout.write(s + '\r')

            # record current time
            self._time = time()

        # Increase time step
        self._timestep += 1

    def Save(self, problem, w, dual=False):
        '''
            Save a variables associated with a time step. Here we assume there
            are two variables where the first variable is vector-valued and the
            second variable is a scalar. If this doesn't fit the particular
            solvers variables the user will need to overload this function.
        '''
        u = w.split()[0]
        p = w.split()[1]

        if self.options['save_frequency'] != 0 \
                and ((self._timestep - 1) % self.options['save_frequency'] == 0
                     or self._timestep == 0 or self._timestep == problem.T):
            if not dual:
                self._ufile << u
                self._pfile << p
            else:
                self._uDualfile << u
                self._pDualfile << p

    def file_naming(self, n=-1, dual=False):
        '''
            Names our files for saving variables. Here we assume there are
            two variables where the first variable is vector-valued and the
            second variable is a scalar. If this doesn't fit the particular
            solvers variables the user will need to overload this function.
        '''
        if n == -1:
            self._ufile = File(self.s + '_u.pvd', 'compressed')
            self._pfile = File(self.s + '_p.pvd', 'compressed')
            self._uDualfile = File(self.s + '_uDual.pvd', 'compressed')
            self._pDualfile = File(self.s + '_pDual.pvd', 'compressed')
            self.meshfile = File(self.s + '_mesh.xml')
        else:
            self._ufile = File(self.s + '_u%02d.pvd' % n, 'compressed')
            self._pfile = File(self.s + '_p%02d.pvd' % n, 'compressed')
            self._uDualfile = File(self.s + '_uDual%02d.pvd' % n, 'compressed')
            self._pDualfile = File(self.s + '_pDual%02d.pvd' % n, 'compressed')
            self.meshfile = File(self.s + '_mesh%02d.xml' % n)

    # this is a separate function so that it can be overloaded
    def Plot(self, problem, W, w):
        '''
            Plots our variables associated with a time step. Here we assume
            there are two variables where the first variable is vector-valued
            and the second variable is a scalar. If this doesn't fit the
            particular solvers variables the user will need to overload this
            function.
        '''
        u = w.split()[0]
        p = w.split()[1]

        if self.vizU is None:
            # Plot velocity and pressure
            self.vizU = plot(u, title='Velocity', rescale=True)
            self.vizP = plot(p, title='Pressure', rescale=True, elevate=0.0)
        else:
            self.vizU.plot(u)
            self.vizP.plot(p)

    def prefix(self, problem):
        '''
            Obtains the beginning of file naming, e.g. Probem Name, Solver Name,
            dimension, etc.
        '''
        # Return file prefix for output files
        p = problem.__module__.split('.')[-1]
        if problem.mesh.topology().dim() > 2:
            p += '3D'
        else:
            p += '2D'
        s = self.__module__.split('.')[-1]
        if(self.options['stabilize']
           and 'stabilization_parameters' in dir(self)):
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

        # Return file suffix for output files
        if not self.options['inviscid'] and self.Re is not None:
            s = 'Re' + str(int(self.Re))
        if self.Ro is not None:
            s += 'Ro' + str(self.Ro)
        if self.Fr is not None:
            s += 'Fr' + str(self.Fr)
        if self.Th is not None:
            s += 'Th' + str(self.Th)

        s += 'T' + str(problem.T) + 'Nx' + str(problem.Nx)
        if problem.Ny is not None:
            s += 'Ny' + str(problem.Ny)
        if problem.mesh.topology().dim() > 2 and problem.Nz is not None:
            s += 'Nz' + str(problem.Nz)

        s += 'K' + str(int(1. / problem.k))

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
