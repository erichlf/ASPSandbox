__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2014-01-26"

from AFES import *
from AFES import Problem as ProblemBase

'''
    This is the implementation for four-gyre problem with double-gyre forcing.
    It is specifically written with the q-psi formulation of QGE in mind.
'''


class Problem(ProblemBase):

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        self.Nx = options['Nx']
        self.Ny = options['Ny']
        self.mesh = RectangleMesh(0, -1, 1, 1, self.Nx, self.Ny)

        # initial time, final time, time-step
        self.t0 = 0.
        self.T = options['T']
        psi_avg = 2.  # expected time-averaged max streamfunction
        self.k = self.time_step(psi_avg, self.mesh)

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

    def initial_conditions(self, W):
        w0 = Expression(('0', '0'))
        w0 = project(w0, W)

        return w0

    def boundary_conditions(self, W, t):
        noslipQ = DirichletBC(W.sub(0), Constant(0.), 'on_boundary')
        noslipPsi = DirichletBC(W.sub(1), Constant(0.), 'on_boundary')

        bcs = [noslipQ, noslipPsi]

        return bcs

    def time_step(self, psi, mesh):
        C_CFL = 100.
        return C_CFL * mesh.hmin()/psi

    def F(self, t):  # Forcing function
        return Expression('sin(pi*x[1])', t=t)

    def functional(self, W, w):  # functional for adaptivity
        (q, psi) = (w[0], w[1])

        M = Constant(0.5) * q * q * dx

        return M

    def __str__(self):
        return 'DoubleGyre'
