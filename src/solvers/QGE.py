__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2013-08-27"

from AFES import *
from AFES import Solver as SolverBase


class Solver(SolverBase):

    '''
        Solver class for solving the Quasi-Geostrophic equations in
        Steamfunction-Vorticity form.
    '''

    def __init__(self, options):
        SolverBase.__init__(self, options)
        try:
            self.stabilize = options['stabilize']
        except:
            self.stabilize = True

    def function_space(self, mesh):
        # Define function spaces
        Q = FunctionSpace(mesh, 'CG', 1)
        P = FunctionSpace(mesh, 'CG', 1)
        W = MixedFunctionSpace([Q, P])

        return W

    def strong_residual(self, w, w2):  # Not implemented yet
        (q, psi) = (w[0], w[1])
        q2 = w2[0]

        R1 = self.Jac(psi, q2) - psi.dx(0)
        R2 = q

        return R1, R2

    def weak_residual(self, problem, k, W, w_theta, w, w_, wt, ei_mode=False):
        (q, psi) = (w_theta[0], w_theta[1])
        Q = w[0]
        Q_ = w_[0]
        (p, chi) = (wt[0], wt[1])

        h = CellSize(W.mesh())  # mesh size

        Re = Constant(problem.Re)  # Reynolds Number
        Ro = Constant(problem.Ro)  # Rossby Number

        t = problem.t0
        f = problem.F(t)  # forcing and mass source/sink

        d1 = conditional(le(h, 1. / Re), h**2, h)  # stabilization parameter
        d2 = 0.01 * conditional(le(h, Ro), h**2, h)  # stabilization parameter

        if(not self.stabilize or ei_mode):
            d1 = Constant(0)
            d2 = Constant(0)

        # weak form of the equations
        r = ((1. / k) * (Q - Q_) * p
             + 1. / Re * inner(grad(q), grad(p))
             + self.Jac(psi, q) * p
             - psi.dx(0) * p) * dx  # vorticity equation
        r -= f * p * dx  # forcing function
        # streamfunction equation
        r += (q * chi - Ro * inner(grad(psi), grad(chi))) * dx

        # least squares stabilization
        # R1, R2 = self.strong_residual(w_theta, w_theta)
        # Rv1, Rv2 = self.strong_residual(wt, w_theta)
        # r += (d1 * (R1 - f) * Rv1 + d2 * R2 * Rv2) * dx
        r += d1 * inner(grad(q), grad(p)) * dx
        r -= d2 * inner(grad(psi), grad(chi)) * dx

        return r

    def Jac(self, psi, q):
        return psi.dx(1) * q.dx(0) - psi.dx(0) * q.dx(1)

    def suffix(self, problem):
        '''
            Obtains the run specific data for file naming, e.g. Nx, k, etc.
        '''

        s = 'Re' + str(problem.Re) + 'Ro' + str(problem.Ro)

        s += 'T' + str(problem.T)
        if problem.Nx is not None:
            s += 'Nx' + str(problem.Nx)
        if problem.Ny is not None:
            s += 'Ny' + str(problem.Ny)

        s += 'K' + str(problem.k)

        return s

    def file_naming(self, problem, n=-1, opt=False):
        s = self.dir + self.prefix(problem) + self.suffix(problem)

        if n == -1:
            if opt:
                self._ufile = File(s + '_qOpt.pvd', 'compressed')
                self._pfile = File(s + '_psiOpt.pvd', 'compressed')
            else:
                self._ufile = File(s + '_q.pvd', 'compressed')
                self._pfile = File(s + '_psi.pvd', 'compressed')
            self._uDualfile = File(s + '_qDual.pvd', 'compressed')
            self._pDualfile = File(s + '_psiDual.pvd', 'compressed')
            self.meshfile = File(s + '_mesh.xml')
        else:
            if self.eifile is None:  # error indicators
                self.eifile = File(s + '_ei.pvd', 'compressed')
            self._ufile = File(s + '_q%d.pvd' % n, 'compressed')
            self._pfile = File(s + '_psi%d.pvd' % n, 'compressed')
            self._uDualfile = File(s + '_qDual%d.pvd' % n, 'compressed')
            self._pDualfile = File(s + '_psiDual%d.pvd' % n, 'compressed')
            self.meshfile = File(s + '_mesh%02d.xml' % n)

    def Save(self, problem, w, dual=False):
        q = w.split()[0]
        psi = w.split()[1]

        if self.saveFrequency != 0 \
                and (self._timestep - 1) % self.saveFrequency == 0:
            if not dual:
                self._ufile << q
                self._pfile << psi
            else:
                self._uDualfile << q
                self._pDualfile << psi

    def Plot(self, problem, W, w, dual=False):
        q = w.split()[0]
        psi = w.split()[1]

        if self.vizU is None:
            self.vizU = plot(q, title='Vorticity', rescale=True, elevate=0.0)
            self.vizP = plot(psi, title='Streamfunction', rescale=True,
                             elevate=0.0)
        else:
            self.vizU.plot(q)
            self.vizP.plot(psi)

    def __str__(self):
        return 'QGE'
