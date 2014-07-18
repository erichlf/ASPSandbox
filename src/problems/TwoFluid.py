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

x0 = 0.
x1 = 1.
y0 = 0.
y1 = 1.

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
        Nx = options['Nx']
        Ny = options['Ny']
        self.mesh = RectangleMesh(x0,y0,x1,y1,Nx, Ny)
        self.Fr = options['Fr']

        self.t0 = 0.
        self.T = options['T']
        self.k = options['dt']

    def initial_conditions(self, V, R, Q):
        U0 = Constant((0, 0))
        rho0 = Expression('x[0]<(x1-x0)/2. ? 1 : 0', x1=x1, x0=x0)
        p0 = Constant(0)

        return U0, rho0, p0

    def boundary_conditions(self, W, t):
        # Create no-slip boundary condition for velocity
        bcs = DirichletBC(W.sub(0), Constant((0.0, 0.0)), NoslipBoundary())

        return bcs

    def F1(self, t):
        #forcing function for the momentum equation
        return Expression(('0.','1./Fr'),Fr=self.Fr,t=t)

    def F2(self, t):
        #mass source for the continuity equation
        return Expression('0.0',t=t)

    def __str__(self):
        return 'TwoFluid'
