__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__ = "GNU GPL version 3 or any later version"

from AFES import *
from AFES import Solver as SolverBase


def f(f0, beta):  # Coriolis parameter
    return Expression('f0 + beta*x[1]', f0=f0, beta=beta)


class Solver(SolverBase):

    '''
        Solver class for solving the Rotational Shallow Water Equation.
    '''

    def __init__(self, options):
        SolverBase.__init__(self, options)
        try:
            self.Re = options['Re']
        except:
            self.Re = 1000
        try:
            self.H = options['H']
        except:
            self.H = 1
        try:
            self.Ro = options['Ro']
        except:
            self.Ro = 1
        try:
            self.Fr = options['Fr']
        except:
            self.Fr = 1
        try:
            self.Th = options['Theta']
        except:
            self.Th = 1
        try:
            linear = options['linear']
        except:
            linear = False
        try:
            inviscid = options['inviscid']
        except:
            inviscid = False

        # if we want a linear version then make a coefficient zero for the
        # terms which only occur in the non-linear from of SWE
        if(linear):
            self.NonLinear = 0
        else:
            self.NonLinear = 1

        if(inviscid):
            self.inviscid = 0
        else:
            self.inviscid = 1

    # strong residual for cG(1)cG(1)
    def strong_residual(self, w, w2):
        (u, eta) = (as_vector((w[0], w[1])), w[2])
        (U, Eta) = (as_vector((w2[0], w2[1])), w2[2])

        # get problem parameters
        Re = self.Re  # Reynolds number
        H = self.H  # Fluid depth
        Ro = self.Ro  # Rossby number
        Th = self.Th  # average wave height
        Fr = self.Fr  # Froude number

        NonLinear = self.NonLinear

        # momentum equation
        R1 = NonLinear * grad(u) * U \
            + 1 / Ro * as_vector((-u[1], u[0])) \
            + Fr ** (-2) * Th * grad(eta)
        # continuity equation
        R2 = 1 / Th * H * div(u)

        return R1, R2

    # weak residual for cG(1)cG(1)
    def weak_residual(self, problem, k, W, w, ww, w_, wt, ei_mode=False):
        (u, eta) = (as_vector((w[0], w[1])), w[2])
        (U, Eta) = (as_vector((ww[0], ww[1])), ww[2])
        (U_, Eta_) = (as_vector((w_[0], w_[1])), w_[2])
        (v, chi) = (as_vector((wt[0], wt[1])), wt[2])

        h = CellSize(W.mesh())  # mesh size
        d1, d2 = self.stabilization_parameters(
            U_, Eta_, k, h)  # stabilization parameters

        # get problem parameters
        Re = self.Re  # Reynolds number
        H = self.H  # Fluid depth
        Ro = self.Ro  # Rossby number
        Th = self.Th  # average wave height
        Fr = self.Fr  # Froude number

        NonLinear = self.NonLinear
        inviscid = self.inviscid

        t0 = problem.t0

        t = t0 + k
        # forcing and mass source/sink
        F1 = problem.F1(t)
        F2 = problem.F2(t)

        # least squares stabilization
        if(ei_mode):
            d1 = 0
            d2 = 0

        # weak form of the equations
        # momentum equation
        r = ((1. / k) * inner(U - U_, v)
                 + 1 / Ro * (u[0] * v[1] - u[1] * v[0])
                 - Fr ** (-2) * Th * eta * div(v)) * dx
        r += inviscid / Re * inner(grad(u), grad(v)) * dx
        # add the terms for the non-linear SWE
        r += NonLinear * inner(grad(u) * u, v) * dx
        # continuity equation
        r += ((1. / k) * (Eta - Eta_) * chi
                  + H / Th * div(u) * chi) * dx

        r -= (inner(F1, v) + F2 * chi) * dx

        R1, R2 = self.strong_residual(w, w)
        Rv1, Rv2 = self.strong_residual(wt, w)
        r += (d1 * inner(R1 - F1, Rv1) + d2 * (R2 - F2) * Rv2) * dx

        return r

    def stabilization_parameters(self, U_, eta_, k, h):
        K1 = (self.Ro * self.Fr ** 2 * self.Th ** (-1)) / 2
        K2 = self.Th / (2 * self.H)
        d1 = K1 * (k ** (-2) + inner(U_, U_) * h ** (-1)) ** (-0.5)
        d2 = K2 * (k ** (-2) + eta_ * eta_ * h ** (-1)) ** (-0.5)

        return d1, d2

    def file_naming(self, n=-1, dual=False):
        if n == -1:
            self._ufile = File(self.s + '_u.pvd', 'compressed')
            self._pfile = File(self.s + '_eta.pvd', 'compressed')
            self._uDualfile = File(self.s + '_uDual.pvd', 'compressed')
            self._pDualfile = File(self.s + '_etaDual.pvd', 'compressed')
            self.meshfile = File(self.s + '_mesh.xml')
        else:
            self._ufile = File(self.s + '_u%d.pvd' % n, 'compressed')
            self._pfile = File(self.s + '_eta%d.pvd' % n, 'compressed')
            self._uDualfile = File(self.s + '_uDual%d.pvd' % n, 'compressed')
            self._pDualfile = File(self.s + '_etaDual%d.pvd' % n, 'compressed')
            self.meshfile = File(self.s + '_mesh%d.xml' % n)

    # this is a separate function so that it can be overloaded
    def Plot(self, problem, W, w):
        u = w.split()[0]
        p = w.split()[1]

        if self.vizU is None:
            # Plot velocity and pressure
            self.vizU = plot(u, title='Velocity', rescale=True)
            self.vizP = plot(p, title='Height', rescale=True)
        else:
            self.vizU.plot(u)
            self.vizP.plot(p)

    def __str__(self):
        return 'SWE'
