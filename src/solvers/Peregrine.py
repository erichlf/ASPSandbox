__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__ = "GNU GPL version 3 or any later version"

from solverbase import *


class Solver(SolverBase):

    '''
        Solver class for Peregrine system.
    '''

    def __init__(self, options):
        SolverBase.__init__(self, options)

        # Parameters
        self.Re = None
        self.Ro = None
        self.Hfile = None
        self.vizH = None

    # define functions spaces
    def function_space(self, mesh):
        '''
            The Peregrine system doesn't quite fit into our normal scheme, so we
            need to define the functions space for the Peregrine System. This is
            mostly so that the wave object is represented.
        '''
        V = VectorFunctionSpace(mesh, 'CG', self.Pu)
        Q = FunctionSpace(mesh, 'CG', self.Pp)
        self.Q = Q  # Bad hack for being able to project into Q space
        W = MixedFunctionSpace([V, Q])

        self.Zeta = Function(Q, name='zeta')
        self.Zeta_ = Function(Q, name='zeta_')
        self.Zeta__ = Function(Q, name='zeta__')
        self.H = Function(Q, name='Bathymetry')
        self.H_ = Function(Q, name='PreviousBathymetry')

        return W

    def strong_residual(self, problem, H, w, wt):
        '''
            Defines the strong residual for Peregrine System.
        '''
        (u, eta) = (as_vector((w[0], w[1])), w[2])
        (v, chi) = (as_vector((wt[0], wt[1])), wt[2])
        # Parameters
        sigma = problem.sigma
        epsilon = problem.epsilon

        k = problem.k

        zeta_tt = 1. / k ** 2 * (self.Zeta - 2 * self.Zeta_ + self.Zeta__)
        zeta_t = 1. / k * (self.Zeta - self.Zeta_)

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
        (u, eta) = (as_vector((w[0], w[1])), w[2])
        (U, Eta) = (as_vector((ww[0], ww[1])), ww[2])
        (U_, Eta_) = (as_vector((w_[0], w_[1])), w_[2])
        (v, chi) = (as_vector((wt[0], wt[1])), wt[2])

        h = CellSize(W.mesh())  # mesh size
        d = 0.1 * h ** (3. / 2.)  # stabilization parameter
        d1, d2 = self.stabilization_parameters(
            U_, Eta_, k, h)  # stabilization parameters

        # set up error indicators
        Z = FunctionSpace(W.mesh(), "DG", 0)
        z = TestFunction(Z)

        theta = self.theta  # time stepping method

        # Parameters
        sigma = problem.sigma
        epsilon = problem.epsilon

        t0 = problem.t0  # initial time
        self.t0 = t0
        T = problem.T  # Final time

        # set up our wave object
        self.wave_object(problem, self.Q, t0, k)

        zeta_tt = 1. / k ** 2 * (self.Zeta - 2 * self.Zeta_ + self.Zeta__)
        zeta_t = 1. / k * (self.Zeta - self.Zeta_)

        # Time stepping method
        zeta_theta = (1. - theta) * self.Zeta_ + theta * self.Zeta
        H_theta = (1. - theta) * self.H_ + theta * self.H

        t = t0 + k

        if(not self.options["stabilize"] or ei_mode):
            d = 0
            d1 = 0
            d2 = 0
        if(not ei_mode):
            z = 1.

        # weak form of the equations
        r = z * (1. / k * inner(U - U_, v) + epsilon * inner(grad(u) * u, v)
                 - div(v) * eta) * dx

        r += z * (sigma ** 2 * 1. / k * div(H_theta * (U - U_)) *
                  div(H_theta * v / 2.) - sigma ** 2 * 1. / k * div(U - U_) *
                  div(H_theta ** 2 * v / 6.)) * dx
        r += z * sigma ** 2 * zeta_tt * div(H_theta * v / 2.) * dx

        r += z * (1. / k * (Eta - Eta_) * chi + zeta_t * chi) * dx
        r -= z * inner(u, grad(chi)) * (epsilon * eta + H_theta) * dx

        r += z * d * \
            (inner(grad(u), grad(v)) + inner(grad(eta), grad(chi))) * dx

        R1, R2, z1, z2 = self.strong_residual(problem, H_theta, w, w)
        Rv1, Rv2, z1, z2 = self.strong_residual(problem, H_theta, w, wt)
        r += z * (d1 * inner(R1 + z1, Rv1) + d2 * (R2 + z2) * Rv2) * dx

        return r

    def stabilization_parameters(self, u, eta, k, h):
        K1 = 2.
        K2 = 2.
        d1 = K1 * (k ** (-2) + inner(u, u) * h ** (-2)) ** (-0.5)
        d2 = K2 * (k ** (-2) + eta * eta * h ** (-2)) ** (-0.5)

        return d1, d2

    def wave_object(self, problem, Q, t, k):
        '''
            Update for wave object at each time step.
        '''
        annotate = self.options['adaptive'] or self.options['optimize']
        H, H_, zeta, zeta_, zeta__ = problem.update_bathymetry(
            Q, t, annotate=annotate)

        self.Zeta.assign(zeta)
        self.Zeta_.assign(zeta_)
        self.Zeta__.assign(zeta__)

        self.H.assign(H)
        self.H_.assign(H_)

    def Save(self, problem, w, dual=False):
        '''
            The Peregrine system doesn't quite fit into our normal scheme, so we
            need to save files differently.
        '''
        u = w.split()[0]
        eta = w.split()[1]

        if self.options['save_frequency'] != 0 and (self._timestep - 1) % self.options['save_frequency'] == 0:
            if not dual:
                self._ufile << u
                self._pfile << eta
                self._Hfile << self.H
            else:
                self._uDualfile << u
                self._pDualfile << eta

    def file_naming(self, n=-1, dual=False):
        '''
            File naming for Peregrine system.
        '''
        if n == -1:
            self._ufile = File(self.s + '_u.pvd', 'compressed')
            self._pfile = File(self.s + '_eta.pvd', 'compressed')
            self._Hfile = File(self.s + '_H.pvd', 'compressed')
            self._uDualfile = File(self.s + '_uDual.pvd', 'compressed')
            self._pDualfile = File(self.s + '_etaDual.pvd', 'compressed')
            self.meshfile = File(self.s + '_mesh.xml')
        else:
            self._ufile = File(self.s + '_u%d.pvd' % n, 'compressed')
            self._pfile = File(self.s + '_eta%d.pvd' % n, 'compressed')
            self._Hfile = File(self.s + '_H%d.pvd' % n, 'compressed')
            self._uDualfile = File(self.s + '_uDual%d.pvd' % n, 'compressed')
            self._pDualfile = File(self.s + '_etaDual%d.pvd' % n, 'compressed')
            self.meshfile = File(self.s + '_mesh%d.xml' % n)

    def Plot(self, problem, W, w):
        '''
            The Peregrine system doesn't quite fit into our normal scheme, so we
            need to plot differently.
        '''
        u = w.split()[0]
        eta = w.split()[1]

        # Plot velocity and height and wave object
        if self.vizU is None:
            self.vizU = plot(u, title='Velocity', rescale=True)
            self.vizP = plot(eta, title='Height', rescale=True)
            self.vizH = plot(self.H, title='Wave Object', rescale=True)
        else:
            self.vizU.plot(u)
            self.vizP.plot(eta)
            self.vizH.plot(self.H)

    def __str__(self):
        return 'Peregrine'
