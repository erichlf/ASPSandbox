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

class InitialConditions(Expression):
    def eval(self,values,x):
        values[0] = 0.
        values[1] = 0.
        if x[0] < 0.5:
            values[2] = 1.
        else:
            values[2] = 0.
        values[3] = 0.
    def value_shape(self):
        return (4,)

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
        self.Nx = options['Nx']
        self.Ny = options['Ny']
        self.mesh = RectangleMesh(x0,y0,x1,y1,self.Nx, self.Ny)
        self.Fr = options['Fr']

        self.t0 = 0.
        self.T = options['T']
        self.k = options['dt']

    #get the initial condition and project it
    def initial_conditions(self, W):
        w0 = InitialConditions()
        w0 = project(w0,W)

        return w0

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
