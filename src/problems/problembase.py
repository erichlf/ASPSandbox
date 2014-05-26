__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
#
#   adapted from problembase.py in nsbench originally developed by
#   Anders Logg <logg@simula.no>
#

from dolfin import *
from numpy import linspace
from math import *

class ProblemBase:
#   Base class for all problems.
    def __init__(self, options):

        # Store options
        self.options = options

        # Parameters must be defined by subclass
        self.mesh    = None

        self.t      = 0
        self.T = options["T"] #final time

        self.bcu    = []
        self.bcp    = []

        self.u0     = None
        self.p0     = None
        self.u      = None
        self.p      = None
        self.solver = None
        self.output_location = ''
        #Parameters for the wave object
        self.h0 = None
        self.hb = None
        self.ad = None
        self.vmax = None
        self.D = None
        self.zeta0 = None
        
    def update_problem(self, t, u, p):
#        Update problem at time t

        # Update state
        self.t = t
        self.u = u
        self.p = p

        # Call problem-specific update
        self.update(t, u, p)

    def update(self, t, u, p):
#        Problem-speficic update at time t
        pass

    def tolerance(self, problem):
#       Return tolerance (used as local convergence criterion).
        if str(problem) == 'Channel':
            return 1e-11
        elif str(problem) == 'Cylinder':
            return 1e-7
        else:
            return 1e-6
