__author__ = "Robin Goix <robin.goix.rg@gmail.com>"
__date__ = "2014-04-22"
__license__  = "GNU GPL version 3 or any later version"
#
#   adapted from Square.py originally developed 
#   by Erich L Foster <erichlf@gmail.com>
#

'''
This is a simple problem with no forcing and on a rectangle.
This is basically a blank problem that we can adapt with optional inputs.
'''

from problembase import *
from numpy import array

# No-slip boundary

class NoslipBoundary(SubDomain):
    def inside(self, x, on_boundary):
	bmarg = 1.e-10 + DOLFIN_EPS
	x0=-5./20.
	x1=10./20.
	y0=-2./20.
	y1=2./20.
        return on_boundary #and \
               #(x[1] < y0 +bmarg or x[1] > y1 - bmarg or \
                #x[0] < x0 + bmarg or x[0] > x1 - bmarg)

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
        lambda0 = float(options["lambda0"])
        self.mesh = RectangleMesh(x0/lambda0, y0/lambda0, x1/lambda0, y1/lambda0, Nx, Ny, 'crossed')

    def initial_conditions(self, V, Q):
        u0 = Expression(("0.0", "0.0"))
        eta0 = Expression("0.0")

        return u0, eta0

    def boundary_conditions(self, V, Q, t):
        # Create no-slip boundary condition for velocity
        bcs = DirichletBC(V, [0.0, 0.0], NoslipBoundary())

        return bcs

    def F1(self, t):
        #forcing function for the momentum equation
        return Expression(self.options['F1'],t=t)

    def F2(self, t):
        #mass source for the continuity equation
        return Expression(self.options['F2'],t=t)

    def __str__(self):
        return 'Square'
