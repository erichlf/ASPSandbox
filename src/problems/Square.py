__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"
#
#   adapted from channel.py in nsbench originally developed 
#   by Kent-Andre Mardal <kent-and@simula.no>
#

'''
This is an extremely boring problem with no forcing and on a square.
This is basically a blank problem that we can adapt with optional inputs.
'''

from problembase import *
from numpy import array

# No-slip boundary
class NoslipBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary 

# Problem definition
class Problem(ProblemBase):
#   2D channel flow.

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Create mesh
        Nx = options["Nx"]
        Ny = options["Ny"]
        x0 = float(options["x0"])
        x1 = float(options["x1"])
        y0 = float(options["y0"])
        y1 = float(options["y1"])
        self.mesh = RectangleMesh(x0,y0,x1,y1,Nx, Ny)

        self.t0 = 0.
        self.T = options['T']
        self.k = options['dt']

    def initial_conditions(self, V, Q):
        u0 = Constant((0, 0))
        eta0 = Constant(0)

        return u0, eta0

    def boundary_conditions(self, V, Q, t):
        # Create no-slip boundary condition for velocity
        bcs = DirichletBC(V, Constant((0.0, 0.0)), NoslipBoundary())

        return bcs

    def F1(self, t):
        #forcing function for the momentum equation
        return Expression(self.options['F1'],t=t)

    def F2(self, t):
        #mass source for the continuity equation
        return Expression(self.options['F2'],t=t)

    def __str__(self):
        return 'Square'
