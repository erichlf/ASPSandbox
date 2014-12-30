__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__ = "GNU GPL version 3 or any later version"

from AFES import *
from AFES import Solver as SolverBase


class Solver(SolverBase):

    '''
        Solver class for Peregrine system.
    '''

    def __init__(self, options):
        SolverBase.__init__(self, options)

        # Parameters
        self.vizH = None

    # define functions spaces
    def function_space(self, mesh):
        '''
            The Peregrine system doesn't quite fit into our normal scheme, so we
            need to define the functions space for the Peregrine System. This is
            mostly so that the wave object is represented.
        '''
        V = VectorFunctionSpace(mesh, 'CG', 1)
        Q = FunctionSpace(mesh, 'CG', 1)
        Nu = FunctionSpace(mesh, 'DG', 1)
        Xi = FunctionSpace(mesh, 'DG', 1)
        W = MixedFunctionSpace([V, Q, Nu, Xi])

        self.w__ = Function(W)

        return W

    def strong_residual(self, problem, w, wt, zeta_t, zeta_tt):
        '''
            Defines the strong residual for Peregrine System.
        '''
        (u, eta, zeta, H) = (as_vector((w[0], w[1])), w[2], w[3], w[4])
        (v, chi, nu, xi) = (as_vector((wt[0], wt[1])), wt[2], wt[3], wt[4])
        # Parameters
        sigma = problem.sigma
        epsilon = problem.epsilon

        # strong form for stabilization
        R1 = epsilon * grad(v) * u + grad(chi)
        z1 = -sigma ** 2 * H / 2. * grad(zeta_tt)

        R2 = div(v * (H + epsilon * eta))
        z2 = zeta_t

        return R1, R2, z1, z2

    def weak_residual(self, problem, k, W, w, ww, w_, wt, ei_mode=False):
        '''
            Defines the weak residual for Peregrine System, including LS
            Stabilization.
        '''
        (u, eta, zeta, H) = (as_vector((w[0], w[1])), w[2], w[3], w[4])
        (U, Eta, Zeta, Height) = (as_vector((ww[0], ww[1])),
                                  ww[2], ww[3], ww[4])
        (U_, Eta_, Zeta_, Height_) = (as_vector((w_[0], w_[1])),
                                      w_[2], w_[3], w_[4])
        (v, chi, nu, xi) = (as_vector((wt[0], wt[1])), wt[2], wt[3], wt[4])

        # initial the n - 2 time-step
        self.w__.assign(w_)
        (U__, Eta__, Zeta__, Height__) = (as_vector((self.w__[0], self.w__[1])),
                                      self.w__[2], self.w__[3], self.w__[4])

        h = CellSize(W.mesh())  # mesh size
        d = 0.1 * h ** (3. / 2.)  # stabilization parameter
        d1, d2 = self.stabilization_parameters(U_, Eta_, k, h)

        # Parameters
        sigma = problem.sigma
        epsilon = problem.epsilon
        beta = problem.beta
        D = problem.D

        zeta_tt = 1. / k ** 2 * (Zeta - 2 * Zeta_ + Zeta__)
        zeta_t = 1. / k * (Zeta - Zeta_)

        if ei_mode:
            d = 0
            d1 = 0
            d2 = 0

        # weak form for peregrine
        r = (1. / k * inner(U - U_, v) + epsilon * inner(grad(u) * u, v)
             - div(v) * eta) * dx

        r += (sigma ** 2 * 1. / k * div(H * (U - U_)) *
              div(H * v / 2.) - sigma ** 2 * 1. / k * div(U - U_) *
              div(H ** 2 * v / 6.)) * dx
        r += sigma ** 2 * zeta_tt * div(H * v / 2.) * dx

        r += (1. / k * (Eta - Eta_) * chi + zeta_t * chi) * dx
        r -= inner(u, grad(chi)) * (epsilon * eta + H) * dx

        # the following will deal with the moving bottom using DG
        alpha = Constant(5.0) # Penalty term

        # Mesh-related functions
        n = FacetNormal(W.mesh())
        h_avg = (h('+') + h('-'))/2

        # ( dot(v, n) + |dot(v, n)| )/2.0
        betan = (dot(beta, n) + abs(dot(beta, n)))/2.0

        r -= dot(beta*zeta, grad(nu)) * dx
        r += dot(betan('+') * zeta('+') - betan('-') * zeta('-'),
                 jump(nu)) * dS + dot(betan * zeta, nu) * ds

        # The following is to describe the bathymetry
        r += (H - (D + epsilon * zeta)) * xi * dx

        # stabilization
        R1, R2, z1, z2 = self.strong_residual(problem, w, w, zeta_t, zeta_tt)
        Rv1, Rv2, z1, z2 = self.strong_residual(problem, w, wt, zeta_t, zeta_tt)
        r += (d1 * inner(R1 + z1, Rv1) + d2 * (R2 + z2) * Rv2) * dx
        # streamline diffusion
        r += d * (inner(grad(u), grad(v)) + inner(grad(eta), grad(chi))) * dx

        return r

    def stabilization_parameters(self, u, eta, k, h):
        K1 = 2.
        K2 = 2.
        d1 = K1 * (k ** (-2) + inner(u, u) * h ** (-2)) ** (-0.5)
        d2 = K2 * (k ** (-2) + eta * eta * h ** (-2)) ** (-0.5)

        return d1, d2

    def post_step(self, problem, t, k, W, w, w_):
        self.w__.assign(w_)

    def Save(self, problem, w, dual=False):
        '''
            The Peregrine system doesn't quite fit into our normal scheme, so we
            need to save files differently.
        '''
        u = w.split()[0]
        eta = w.split()[1]
        zeta = w.split()[2]
        H = w.split()[3]

        if self.saveFrequency != 0 and (self._timestep - 1) \
                % self.saveFrequency == 0:
            if not dual:
                u.rename("Velocity", "Peregrine")
                eta.rename("WaveHeight", "Peregrine")
                u.rename("Bathymetry", "Peregrine")
                self._ufile << u
                self._pfile << eta
                self._Hfile << H
            else:
                u.rename("DualVelocity", "Peregrine")
                eta.rename("DualWaveHeight", "Peregrine")
                u.rename("DualBathymetry", "Peregrine")
                self._uDualfile << u
                self._pDualfile << eta
                self._HDualfile << H

    def file_naming(self, problem, n=-1, opt=False):
        '''
            File naming for Peregrine system.
        '''
        s = 'results/' + self.prefix(problem) + self.suffix(problem)

        if n == -1:
            if opt:
                self._ufile = File(s + '_uOpt.pvd', 'compressed')
                self._pfile = File(s + '_etaOpt.pvd', 'compressed')
                self._Hfile = File(s + '_HOpt.pvd', 'compressed')
            else:
                self._ufile = File(s + '_u.pvd', 'compressed')
                self._pfile = File(s + '_eta.pvd', 'compressed')
                self._Hfile = File(s + '_H.pvd', 'compressed')
            self._uDualfile = File(s + '_uDual.pvd', 'compressed')
            self._pDualfile = File(s + '_etaDual.pvd', 'compressed')
            self._pDualfile = File(s + '_HDual.pvd', 'compressed')
            self.meshfile = File(s + '_mesh.xml')
        else:
            self._ufile = File(s + '_u%d.pvd' % n, 'compressed')
            self._pfile = File(s + '_eta%d.pvd' % n, 'compressed')
            self._Hfile = File(s + '_H%d.pvd' % n, 'compressed')
            self._uDualfile = File(s + '_uDual%d.pvd' % n, 'compressed')
            self._pDualfile = File(s + '_etaDual%d.pvd' % n, 'compressed')
            self._pDualfile = File(s + '_HDual%d.pvd' % n, 'compressed')
            self.meshfile = File(s + '_mesh%d.xml' % n)

    def Plot(self, problem, W, w):
        '''
            The Peregrine system doesn't quite fit into our normal scheme, so we
            need to plot differently.
        '''
        u = w.split()[0]
        eta = w.split()[1]
        zeta = w.split()[2]
        H = w.split()[3]

        # Plot velocity and height and wave object
        if self.vizU is None:
            self.vizU = plot(u, title='Velocity', rescale=True)
            self.vizP = plot(eta, title='Height', rescale=True)
            self.vizH = plot(H, title='Bathymetry', rescale=True)
        else:
            self.vizU.plot(u)
            self.vizP.plot(eta)
            self.vizH.plot(H)

    def __str__(self):
        return 'Peregrine'
