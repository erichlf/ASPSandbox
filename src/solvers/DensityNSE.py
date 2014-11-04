__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__ = "GNU GPL version 3 or any later version"

from AFES import *
from AFES import Solver as SolverBase


class Solver(SolverBase):

    '''
        Solver class for solving the variable density incompressible
        Navier-Stokes equation.
    '''

    def __init__(self, options):

        if options['dim'] is not 2:
            options['dim'] = 2

        try:
            self.nu = 1 / options['Re']
        except:
            self.nu = 1E-3
            options['Re'] = 1000

        SolverBase.__init__(self, options)

        self.vizRho = None
        self.rhofile = None

    # Define function spaces
    def function_space(self, mesh):
        V = VectorFunctionSpace(mesh, 'CG', self.Pu)
        R = FunctionSpace(mesh, 'CG', self.Pp)
        Q = FunctionSpace(mesh, 'CG', self.Pp)

        W = MixedFunctionSpace([V, R, Q])

        return W

    # strong residual for cG(1)cG(1)
    def strong_residual(self, w, w2):  # U, v, rho, r, p):
        (u, rho, p) = (as_vector((w[0], w[1])), w[2], w[3])
        (U, Rho, P) = (as_vector((w2[0], w2[1])), w2[2], w2[3])

        R1 = Rho * grad(U) * u + grad(p)
        R2 = div(Rho * u)
        R3 = div(u)

        return R1, R2, R3

    # weak residual for cG(1)cG(1)
    def weak_residual(self, problem, k, W, w, ww, w_, wt, ei_mode=False):
        # rho = 1/rho' - 1
        (u, rho, p) = (as_vector((w[0], w[1])), w[2], w[3])
        (U, Rho, P) = (as_vector((ww[0], ww[1])), ww[2], ww[3])
        (U_, Rho_, P_) = (as_vector((w_[0], w_[1])), w_[2], w_[3])
        (v, psi, q) = (as_vector((wt[0], wt[1])), wt[2], wt[3])

        h = CellSize(W.mesh())  # mesh size
        # stabilization parameters
        d1, d2, d3 = self.stabilization_parameters(U_, Rho_, P_, k, h)

        nu = self.nu  # kinematic viscosity

        t0 = problem.t0

        t = t0 + k
        # forcing and mass source/sink
        f = problem.F(t)

        # least squares stabilization
        if ei_mode:
            d1 = 0
            d2 = 0
            d3 = 0

        # weak form of the equations
        r = ((1. / k) * (Rho - Rho_) + div(rho * u)) * \
            psi * dx  # mass equation
        r += (rho * ((1. / k) * inner(U - U_, v)
                     + inner(grad(u) * u, v))
              + inner(grad(p), v)
              + nu * inner(grad(u), grad(v))
              - rho * dot(f, v)) * dx  # momentum equation
        r += div(u) * q * dx  # continuity equation

        r += d1 * (inner(grad(u), grad(v))) * dx
        r += d2 * (inner(grad(rho), grad(psi))) * dx
        r += d3 * ((inner(grad(p), grad(q))) + p * q) * dx
        return r

    def stabilization_parameters(self, u, rho, p, k, h):
        K1 = 1000
        K2 = 1
        K3 = 1E-4
        d1 = K1 * h ** (2.)
        d2 = K2 * h ** (2.)
        d3 = K3 * h ** (2.)

        return d1, d2, d3

    def Save(self, problem, w, dual=False):
        u = w.split()[0]
        rho = w.split()[1]
        p = w.split()[2]

        if self.options['save_frequency'] != 0 \
                and (self._timestep - 1) % self.options['save_frequency'] == 0:
            if not dual:
                self._ufile << u
                self._pfile << p
                self._rhofile << rho
            else:
                self._uDualfile << u
                self._pDualfile << p
                self._rhoDualfile << rho

    def file_naming(self, n=-1, dual=False):
        if n == -1:
            self._ufile = File(self.s + '_u.pvd', 'compressed')
            self._rhofile = File(self.s + '_rho.pvd', 'compressed')
            self._pfile = File(self.s + '_p.pvd', 'compressed')
            self._uDualfile = File(self.s + '_uDual.pvd', 'compressed')
            self._rhoDualfile = File(self.s + '_rhoDual.pvd', 'compressed')
            self._pDualfile = File(self.s + '_pDual.pvd', 'compressed')
            self.meshfile = File(self.s + '_mesh.xml')
        else:
            self._ufile = File(self.s + '_u%d.pvd' % n, 'compressed')
            self._rhofile = File(self.s + '_rho%d.pvd' % n, 'compressed')
            self._pfile = File(self.s + '_p%d.pvd' % n, 'compressed')
            self._uDualfile = File(self.s + '_uDual%d.pvd' % n, 'compressed')
            self._rhoDualfile = File(
                self.s +
                '_rhoDual%d.pvd' %
                n,
                'compressed')
            self._pDualfile = File(self.s + '_pDual%d.pvd' % n, 'compressed')
            self.meshfile = File(self.s + '_mesh%02d.xml' % n)

    def Plot(self, problem, W, w):
        u = w.split()[0]
        rho = w.split()[1]
        p = w.split()[2]

        if self.vizU is None:
            # Plot velocity and pressure
            self.vizU = plot(u, title='Velocity', rescale=True)
            self.vizP = plot(p, title='Pressure', rescale=True, elevate=0.0)
            self.vizRho = plot(rho, title='Density', rescale=True, elevate=0.0)
        else:
            self.vizU.plot(u)
            self.vizP.plot(p)
            self.vizRho.plot(rho)

    def __str__(self):
        return 'DensityNSE'
