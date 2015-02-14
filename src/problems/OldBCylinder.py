__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2014-08-12"
__copyright__ = "Copyright (C) 2013-2014 " + __author__
__license__ = "GNU GPL version 3 or any later version"
#
#   adapted from ns in nsbench originally developed by
#   Anders Logg <logg@simula.no>
#

from Cylinder import *
from Cylinder import Problem as Cylinder
import sys


class InitialConditions2D(Expression):

    def eval(self, values, x):
        values[0] = 0.
        values[1] = 0.
        values[2] = 0.
        values[3] = 0.
        values[4] = 0.
        values[5] = 0.
        values[6] = 0.

    def value_shape(self):
        return (7,)


class InitialConditions3D(Expression):

    def eval(self, values, x):
        values[0] = 0.
        values[1] = 0.
        values[2] = 0.
        values[3] = 0.
        values[4] = 0.
        values[5] = 0.
        values[6] = 0.
        values[7] = 0.

    def value_shape(self):
        return (8,)


class Problem(Cylinder):

    '''
        Purely a container to call Cylinder. This
        will tell Cylinder.py to use a viscoelastic problem.
    '''

    def __init__(self, options, cube=False):
        Cylinder.__init__(self, options, cube=False)

    def initial_conditions(self, W):
        if self.dim == 2:
            w0 = InitialConditions2D()
        else:
            w0 = InitialConditions3D()

        w0 = project(w0, W)

        return w0

    def boundary_conditions(self, W, t):
        self.U.t = t
        # Create inflow boundary condition
        bc0 = DirichletBC(W.sub(0), self.U, InflowBoundary())

        # Create no-slip boundary condition
        bc1 = DirichletBC(W.sub(0), self.noSlip,
                          NoSlipBoundary(self.dim))

        # Create outflow boundary condition for pressure
        # bc2 = DirichletBC(W.sub(1), Constant(0), OutflowBoundary())

        # Collect boundary conditions
        bcs = [bc0, bc1]

        return bcs

    def F(self, t):
        # forcing function for the momentum equation
        if self.dim == 2:
            f = Constant((0, 0))
        else:
            f = Constant((0, 0, 0))

        return f

    def functional(self, W, w):
        '''
            This is the functional used for adaptivity.
            We assume the problem is much like NSE.
        '''
        if W.mesh().topology().dim() == 2:
            (u, p, tau) = (as_vector((w[0], w[1])), w[2],
                        as_tensor(((w[3], w[4]), (w[5], w[6]))))
        else:
            print
            print 'This problem can only be used in 2D.'
            sys.exit(1)

        # n = FacetNormal(mesh)
        # marker = problem.marker()

        M = u[0] * dx  # Mean of the x-velocity in the whole domain
        # M = marker*p*n[0]*ds  # Drag (only pressure)

        return M

    def __str__(self):
        return 'OldBCylinder'
