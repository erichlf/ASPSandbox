__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
#
#   adapted from problembase.py in nsbench originally developed by
#   Anders Logg <logg@simula.no>
#

from dolfin import *
from dolfin_adjoint import *
from numpy import linspace
from math import *


class Problem:
#   Base class for all problems.

    def __init__(self, options):

        # Store options
        self.options = options

        # Parameters must be defined by subclass
        self.mesh = None

        self.t = 0
        self.T = options["T"]  # final time

        # reset our discretization
        self.Nx = None
        self.Ny = None
        self.Nz = None

        self.solver = None
        self.output_location = ''
