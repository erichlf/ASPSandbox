from AFES import *
from AFES import Solver as SolverBase


class Solver(SolverBase):
    def __init__(self, options):
        SolverBase.__init__(self, options)

        self.eta = Constant(1.)

    # Define function spaces
    def function_space(self, mesh):
        V = VectorFunctionSpace(mesh, 'CG', 1)
        Q = FunctionSpace(mesh, 'CG', 1)
        S = TensorFunctionSpace(mesh, 'CG', 1)

        W = MixedFunctionSpace([V, Q, S])

        return W

    def strong_residual(self, w, w2):
        (u, p, tau) = (as_vector((w[0], w[1])), w[2],
                       as_tensor(((w[3], w[4]), (w[5], w[6]))))

        R1 = tau - 2. * self.eta * self.epsilon(u)
        R2 = - div(tau) + grad(p)
        R3 = div(u)

        return R1, R2, R3

    def weak_residual(self, problem, k, W, w_theta, w, w_, wt, ei_mode=False):

        (u, p, tau) = (as_vector((w_theta[0], w_theta[1])), w_theta[2],
                       as_tensor(((w_theta[3], w_theta[4]),
                                  (w_theta[5], w_theta[6]))))
        (U, T) = (as_vector((w[0], w[1])),
                  as_tensor(((w[3], w[4]), (w[5], w[6]))))
        (U_, T_) = (as_vector((w_[0], w_[1])),
                    as_tensor(((w_[3], w_[4]), (w_[5], w_[6]))))
        (v, q, s) = (as_vector((wt[0], wt[1])), wt[2],
                     as_tensor(((wt[3], wt[4]), (wt[5], wt[6]))))

        eta = self.eta
        nu = problem.nu

        h = CellSize(W.mesh())
        d1 = conditional(le(h, nu), h**2, h)  # stabilization parameter

        if ei_mode:  # turn off stabilization in ei_mode
            d1 = Constant(0)

        f = problem.F(problem.t0)  # forcing

        # Weak form
        r = 1. / k * inner(T - T_, s) * dx
        r += 1. / (2. * eta) * inner(tau - self.epsilon(u), s) * dx
        r += 1. / k * inner(U - U_, v) * dx
        r += (inner(tau, self.epsilon(v)) - p*div(v) - inner(f, v)) * dx
        r += div(u) * q * dx

        # GLS stabilization
        R1, R2, R3 = self.strong_residual(w_theta, w_theta)
        Rv1, Rv2, Rv3 = self.strong_residual(wt, w_theta)
        r += d1 * (inner(R1, Rv1) + inner(R2 - f, Rv2) + R3 * Rv3) * dx

        return r

    def epsilon(self, u):
        return Constant(0.5) * (nabla_grad(u) + nabla_grad(u).T)

    def suffix(self, problem):
        '''
            Obtains the run specific data for file naming, e.g. Nx, etc.
        '''

        try:
            s = 'Re' + str(problem.Re)
        except:
            s = 'nu' + str(problem.nu)

        if problem.Nx is not None:
            s += 'Nx' + str(problem.Nx)
        if problem.Ny is not None:
            s += 'Ny' + str(problem.Ny)

        return s

    def file_naming(self, problem, n=-1, opt=False):
        s = self.dir + self.prefix(problem) + self.suffix(problem)

        if n == -1:
            if opt:
                self._ufile = File(s + '_uOpt.pvd', 'compressed')
                self._taufile = File(s + '_tauOpt.pvd', 'compressed')
                self._pfile = File(s + '_pOpt.pvd', 'compressed')
            else:
                self._ufile = File(s + '_u.pvd', 'compressed')
                self._taufile = File(s + '_tau.pvd', 'compressed')
                self._pfile = File(s + '_p.pvd', 'compressed')
            self._uDualfile = File(s + '_uDual.pvd', 'compressed')
            self._tauDualfile = File(s + '_tauDual.pvd', 'compressed')
            self._pDualfile = File(s + '_pDual.pvd', 'compressed')
            self.meshfile = File(s + '_mesh.xml')
        else:
            if self.eifile is None:  # error indicators
                self.eifile = File(s + '_ei.pvd', 'compressed')
            self._ufile = File(s + '_u%d.pvd' % n, 'compressed')
            self._taufile = File(s + '_tau%d.pvd' % n, 'compressed')
            self._pfile = File(s + '_p%d.pvd' % n, 'compressed')
            self._uDualfile = File(s + '_uDual%d.pvd' % n, 'compressed')
            self._tauDualfile = File(s + '_tauDual%d.pvd' % n, 'compressed')
            self._pDualfile = File(s + '_pDual%d.pvd' % n, 'compressed')
            self.meshfile = File(s + '_mesh%02d.xml' % n)

    def Save(self, problem, w, dual=False):
        u = w.split()[0]
        p = w.split()[1]
        tau = w.split()[2]

        if self.saveFrequency != 0 \
                and (self._timestep - 1) % self.saveFrequency == 0:
            if not dual:
                self._ufile << u
                self._pfile << p
                self._taufile << tau
            else:
                self._uDualfile << u
                self._pDualfile << p
                self._tauDualfile << tau

    def Plot(self, problem, W, w):
        u = w.split()[0]
        p = w.split()[1]
        # tau = w.split()[2]

        if self.vizU is None:
            # Plot velocity and pressure
            self.vizU = plot(u, title='Velocity', rescale=True)
            self.vizP = plot(p, title='Pressure', rescale=True, elevate=0.0)
            # self.vizTau = plot(tau, title='Stress', rescale=True, elevate=0.0)
        else:
            self.vizU.plot(u)
            self.vizP.plot(p)
            # self.vizTau.plot(tau)

    def __str__(self):
        return 'ThreeField_Stokes'
