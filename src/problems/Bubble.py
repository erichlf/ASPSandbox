__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__ = "GNU GPL version 3 or any later version"

from TwoFluid import *
from TwoFluid import Problem as TwoFluid

xcenter, ycenter = 0, d
r = 0.25 * d


class InitialConditions(Expression):

    def __init__(self, At):
        self.rhoMin = 1.
        self.rhoMax = self.rhoMin * (1. + At) / (1. - At)

    def eval(self, values, x):
        rhoMin, rhoMax = self.rhoMin, self.rhoMax

        radius = sqrt((x[0] - xcenter)**2. + (x[1] - ycenter)**2)

        values[0], values[1], values[3] = 0., 0., 0.
        if radius < r + DOLFIN_EPS:
            values[2] = rhoMax
        else:
            values[2] = rhoMin

    def value_shape(self):
        return (4,)


class Problem(TwoFluid):  # This problem is really just the TwoFluid problem

    def __init__(self, options):

        TwoFluid.__init__(self, options)

    def initial_conditions(self, W):

        # artificial viscosity for stabilization
        self.artificial_viscosity(W)

        w0 = InitialConditions(self.At)
        w0 = project(w0, W)

        return w0

    def __str__(self):
        return 'Bubble'
