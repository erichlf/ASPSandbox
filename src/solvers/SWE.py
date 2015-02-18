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
    def strong_residual(self, w, w2, Ro, Fr, Th, H):
        (u, eta) = (as_vector((w[0], w[1])), w[2])
        u2 = as_vector((w2[0], w2[1]))

        NonLinear = self.NonLinear

        # momentum equation
        R1 = NonLinear * grad(u) * u2 \
            + 1 / Ro * as_vector((-u[1], u[0])) \
            + Fr ** (-2) * Th * grad(eta)
        # continuity equation
        R2 = 1 / Th * H * div(u)

        return R1, R2

    # weak residual for cG(1)cG(1)
    def weak_residual(self, problem, k, W, w_theta, w, w_, wt, ei_mode=False):
        (u, eta) = (as_vector((w_theta[0], w_theta[1])), w_theta[2])
        (U, Eta) = (as_vector((w[0], w[1])), w[2])
        (U_, Eta_) = (as_vector((w_[0], w_[1])), w_[2])
        (v, chi) = (as_vector((wt[0], wt[1])), wt[2])

        h = CellSize(W.mesh())  # mesh size
        d1, d2 = self.stabilization_parameters(U_, Eta_, k, h)

        Re = problem.Re
        H = problem.H
        Ro = problem.Ro
        Fr = porblem.Fr
        Th = problem.Theta

        NonLinear = self.NonLinear
        inviscid = self.inviscid

        t = problem.t0

        # forcing and mass source/sink
        F1 = problem.F1(t)
        F2 = problem.F2(t)

        # least squares stabilization
        if(ei_mode):
            d1 = Constant(0)
            d2 = Constant(0)

        # weak form of the equations
        # momentum equation
        r = ((1. / k) * inner(U - U_, v)
             + 1 / Ro * (u[0] * v[1] - u[1] * v[0])
             - Fr ** (-2) * Th * eta * div(v)) * dx
        r += inviscid / Re * inner(grad(u), grad(v)) * dx
        # add the terms for the non-linear SWE
        r += NonLinear * inner(grad(u) * u, v) * dx
        # continuity equation
        r += ((1. / k) * (Eta - Eta_) * chi + H / Th * div(u) * chi) * dx

        r -= (inner(F1, v) + F2 * chi) * dx

        R1, R2 = self.strong_residual(w_theta, w_theta, Ro, Fr, Th, H)
        Rv1, Rv2 = self.strong_residual(wt, w_theta, Ro, Fr, Th, H)
        r += (d1 * inner(R1 - F1, Rv1) + d2 * (R2 - F2) * Rv2) * dx

        return r

    def stabilization_parameters(self, U_, eta_, k, h):
        K1 = (self.Ro * self.Fr ** 2 * self.Th ** (-1)) / 2
        K2 = self.Th / (2 * self.H)
        d1 = K1 * (k ** (-2) + inner(U_, U_) * h ** (-1)) ** (-0.5)
        d2 = K2 * (k ** (-2) + eta_ * eta_ * h ** (-1)) ** (-0.5)

        return d1, d2

    def suffix(self, problem):
        '''
            Obtains the run specific data for file naming, e.g. Nx, k, etc.
        '''

        s = 'Re' + str(problem.Re) + 'Fr' + str(problem.Fr) \
            + 'H' + str(problem.H) + 'Theta' + str(problem.Th)

        s += 'T' + str(problem.T)
        if problem.Nx is not None:
            s += 'Nx' + str(problem.Nx)
        if problem.Ny is not None:
            s += 'Ny' + str(problem.Ny)
        if problem.mesh.topology().dim() > 2 and problem.Nz is not None:
            s += 'Nz' + str(problem.Nz)

        s += 'K' + str(problem.k)

        return s

    def file_naming(self, problem, n=-1, opt=False):

        s = self.dir + self.prefix(problem) + self.suffix(problem)

        if n == -1:
            if opt:
                self._ufile = File(s + '_uOpt.pvd', 'compressed')
                self._pfile = File(s + '_etaOpt.pvd', 'compressed')
            else:
                self._ufile = File(s + '_u.pvd', 'compressed')
                self._pfile = File(s + '_eta.pvd', 'compressed')
            self._uDualfile = File(s + '_uDual.pvd', 'compressed')
            self._pDualfile = File(s + '_etaDual.pvd', 'compressed')
            self.meshfile = File(s + '_mesh.xml')
        else:
            self._ufile = File(s + '_u%d.pvd' % n, 'compressed')
            self._pfile = File(s + '_eta%d.pvd' % n, 'compressed')
            self._uDualfile = File(s + '_uDual%d.pvd' % n, 'compressed')
            self._pDualfile = File(s + '_etaDual%d.pvd' % n, 'compressed')
            self.meshfile = File(s + '_mesh%d.xml' % n)

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
