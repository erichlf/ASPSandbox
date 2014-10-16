__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
#
#   adapted from problembase.py in nsbench originally developed by
#   Anders Logg <logg@simula.no>
#

from dolfin import *
from dolfin_adjoint import *
from math import *


class ProblemBase:  # Base class for all problems.

    def __init__(self, options):

        # Store options
        self.options = options

        # Parameters must be defined by subclass
        self.mesh = None

        # time domain and time step
        self.t0 = 0
        self.t = self.t0
        self.T = options["T"]  # final time
        self.k = options['k']

        # reset our discretization
        self.Nx = None
        self.Ny = None
        self.Nz = None

        self.solver = None
        self.output_location = ''

    def update(self, W, t):
        pass
