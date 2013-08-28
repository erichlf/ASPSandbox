__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2013-08-27"
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

    def solve(self, problem, dt, plot_solution=True):
#       Solve problem
        raise NotImplementedError

    def prefix(self, problem):
        #Return file prefix for output files
        p = problem.__module__.split('.')[-1].lower()
        s = self.__module__.split('.')[-1].lower()
        return problem.output_location + p + "_" + s

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
        if self.options["save_solution"]:
            # Save velocity and pressure
            frequency = self.options["save_frequency"]
            N = self.options["N"]
            if (self._timestep - 1) % frequency == 0:
                # Create files for saving
                if self._ufile is None:
                    self._ufile = File('results/' + self.prefix(problem) +'N'+ \
                            str(N) + '_u.pvd')
                if self._pfile is None:
                    self._pfile = File('results/' + self.prefix(problem) +'N'+ \
                            str(N) + '_p.pvd')
                self._ufile << u
                self._pfile << p

        # Plot solution
        if self.options["plot_solution"]: 
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
                self.vizP.plot(p)

        # Check memory usage
        if self.options["check_mem_usage"]:
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
