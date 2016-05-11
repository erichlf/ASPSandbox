__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2014-01-26"

from ASP import *
from ASP import Problem as ProblemBase

'''
    This is the implementation for four-gyre problem with double-gyre forcing.
    It is specifically written with the q-psi formulation of QGE in mind.
'''


class Problem(ProblemBase):

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        try:
            self.potential = options['potential_vorticity']
        except:
            self.potential = False

        self.Nx = options['Nx']
        self.Ny = options['Ny']
        self.mesh = RectangleMesh(Point(0, -1), Point(1, 1), self.Nx, self.Ny)

        # initial time, final time, time-step
        self.t0 = 0.
        self.T = options['T']
        self.Ubar = 2.  # expected time-averaged max streamfunction
        try:
            self.k = options['k']
        except:
            self.k = 0.01
        # self.k = self.time_step(self.Ubar, self.mesh)

        # Reynolds number
        try:
            self.Re = options['Re']
        except:
            self.Re = 200

        # Rossby number
        try:
            self.Ro = options['Ro']
        except:
            self.Ro = 0.0016

    def initial_conditions(self, W, annotate=False):
        if(self.potential):
            w0 = Expression(('1. / Ro * x[1]', '0'), Ro=self.Ro)
        else:
            w0 = Expression(('0', '0'))
        w0 = project(w0, W)

        return w0

    def boundary_conditions(self, W, t):
        if(self.potential):
            qBC = Expression('1. / Ro * x[1]', Ro=self.Ro)
        else:
            w0 = Expression(('0', '0'))
        noslipQ = DirichletBC(W.sub(0), qBC, 'on_boundary')
        noslipPsi = DirichletBC(W.sub(1), Constant(0.), 'on_boundary')

        bcs = [noslipQ, noslipPsi]

        return bcs

    # def time_step(self, psi, mesh):
    #     C_CFL = 10.
    #     return 0.01  # C_CFL * mesh.hmin()/psi

    def F(self, t):  # Forcing function
        return Expression('sin(pi*x[1])', t=t)

    def functional(self, W, w):  # functional for adaptivity
        q = w[0]

        M = Constant(0.5) * q * q * dx

        return M

    def __str__(self):
        return 'DoubleGyre'
