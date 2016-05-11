__author__ = "Erich L Foster <efoster@bcamath.org>"
__date__ = "2016-06-11"

from ASP import *
from ASP import Solver as SolverBase


class Solver(SolverBase):

    '''
        Solver class for solving the Quasi-Geostrophic equations in
        Steamfunction-Vorticity form.
    '''

    def __init__(self, options):
        SolverBase.__init__(self, options)
        try:
            self.dg = options['dg']
        except:
            self.dg = True

    def function_space(self, mesh):
        # Define function spaces
        if(self.dg):
            Q = FunctionSpace(mesh, 'DG', 1)
        else:
            Q = FunctionSpace(mesh, 'CG', 1)
        P = FunctionSpace(mesh, 'CG', 1)
        W = MixedFunctionSpace([Q, P])

        return W

    def weak_residual(self, problem, k, W, w_theta, w, w_, wt, ei_mode=False):
        (q, psi) = (w_theta[0], w_theta[1])
        Q = w[0]
        Q_ = w_[0]
        (p, chi) = (wt[0], wt[1])

        b = as_vector((psi.dx(1), - psi.dx(0)))

        y = Expression('x[1]')  # y-coordinate

        Re = Constant(problem.Re)  # Reynolds Number
        Ro = Constant(problem.Ro)  # Rossby Number

        t = problem.t0
        f = problem.F(t)  # forcing and mass source/sink

        if(self.dg):
            flux = self.Flux(problem, W, q, psi, p, chi)
        else:
            flux = 0 * dx

        # weak form of the equations
        r = ((1. / k) * (Q - Q_) * p
             + 1. / Re * inner(grad(q), grad(p))
             + dot(b, grad(q)) * p) * dx # vorticity equation
        r -= f * p * dx  # forcing function
        r += flux # flux for dg method

        # streamfunction equation
        r += (q * chi + Ro * inner(grad(psi), grad(chi)) - y * chi) * dx

        return r

    def Flux(self, problem, W, q, psi, p, chi):
        Re = Constant(problem.Re)
        n = FacetNormal(W.mesh())
        h = CellSize(W.mesh())  # mesh size

        b = as_vector((psi.dx(1), - psi.dx(0)))
        bn = (dot(b, n) + abs(dot(b, n))) / 2.

        alpha = Constant(5.)

        flux_fac = 1. / Re * (alpha / h('+')) * dot(jump(p, n), jump(q, n)) * dS \
                   - 1. / Re * dot(avg(grad(p)), jump(q, n)) * dS \
                   - 1. / Re * dot(jump(p, n), avg(grad(q))) * dS

        flux_vel = dot(jump(p), bn('+') * q('+') - bn('-') * q('-')) * dS \
                   + dot(p, bn * q) * ds

        return flux_fac + flux_vel

    def suffix(self, problem):
        '''
            Obtains the run specific data for file naming, e.g. Nx, k, etc.
        '''

        if(self.dg):
            s = '-DG-'
        else:
            s = '-CG-'

        s += 'Re' + str(problem.Re) + 'Ro' + str(problem.Ro)

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
        return 'QGECGDG'
