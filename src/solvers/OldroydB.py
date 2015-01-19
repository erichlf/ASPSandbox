__author__ = "Erich L Foster <erichlf@gmail.com>, Susanne Claus <susanne.claus@ucl.ac.uk>"
__date__ = "2015-01-14"
__license__ = "GNU GPL version 3 or any later version"

from AFES import *
from AFES import Solver as SolverBase


class Solver(SolverBase):

    '''
        Solver class for solving viscoelastic Oldroyd-B equation.
    '''

    def __init__(self, options):

        SolverBase.__init__(self, options)

        self.vizRho = None
        self.rhofile = None

    # Define function spaces
    def function_space(self, mesh):
        V = VectorFunctionSpace(mesh, 'CG', 1)
        Q = FunctionSpace(mesh, 'CG', 1)
        S = TensorFunctionSpace(mesh, 'CG', 1)

        W = MixedFunctionSpace([V, Q, S])

        return W

    # strong residual for cG(1)cG(1)
    def strong_residual(self, w, w2):  # U, v, rho, r, p):
        (u, p, tau) = (as_vector((w[0], w[1])), w[2], as_tensor((w[3], w[4]),(w[5], w[6])))
        (U, P, T) = (as_vector((w2[0], w2[1])), w2[2], as_tensor((w2[3], w2[4]),(w2[5], w2[6])))

        R1 = grad(U) * u + grad(p)
        R2 = div( u)
        R3 = div(u)

        return R1, R2, R3

    def epsilon(self,u):
        return 0.5*(nabla_grad(u) + nabla_grad(u).T)

    def a_h(self,sigma, v):
        return inner(sigma,self.epsilon(v))*dx

    def b_h(self,p,v):
        return -div(v)*p*dx

    # weak residual for cG(1)cG(1)
    def weak_residual(self, problem, k, W, w, ww, w_, wt, ei_mode=False):
        (u, p, tau) = (as_vector((w[0], w[1])), w[2], as_tensor(((w[3], w[4]),(w[5], w[6]))))
        (U, P, T) = (as_vector((ww[0], ww[1])), ww[2], as_tensor(((ww[3], ww[4]),(ww[5], ww[6]))))
        (U_, P_, T_) = (as_vector((w_[0], w_[1])), w_[2], as_tensor(((w_[3], w_[4]),(w_[5], w_[6]))))
        (v, q, s) = (as_vector((wt[0], wt[1])), wt[2], as_tensor(((wt[3], wt[4]),(wt[5], wt[6]))))

        h = CellSize(W.mesh())  # mesh size
        # stabilization parameters
        d1, d2 = self.stabilization_parameters(U_, P_, T_, k, h)

        nu = problem.nu  # kinematic viscosity
        rho = Constant(1)  # kinematic viscosity

        t0 = problem.t0

        t = t0 + k
        # forcing and mass source/sink
        f = problem.F1(t)

        # least squares stabilization
        if ei_mode:
            d1 = 0
            d2 = 0
            d3 = 0

        # weak form of the equations
        r = rho * ((1. / k) * inner(U - U_, v)
                     + inner(grad(u) * u, v))*dx
        r += 1/(2*nu)*inner(tau,s)*dx
        r += self.a_h(tau, v) - self.a_h(s, u)
        r += self.b_h(p, v) - self.b_h(q,u)
        r +=  - rho * dot(f, v)* dx  # momentum equation
        r += div(u) * q * dx  # continuity equation

        r += d1 * (inner(grad(u), grad(v))) * dx
        r += d2 * ((inner(grad(p), grad(q))) + p * q) * dx
        return r

    def stabilization_parameters(self, u, p, tau, k, h):
        K1 = 1000
        K2 = 1E-1
        d1 = K1 * h ** (2.)
        d2 = K2 * h ** (2.)

        return d1, d2

    def suffix(self, problem):
        '''
            Obtains the run specific data for file naming, e.g. Nx, k, etc.
        '''

        s = 'nu' + str(problem.nu)

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
        s = 'results/' + self.prefix(problem) + self.suffix(problem)

        if n == -1:
            if opt:
                self._ufile = File(s + '_uOpt.pvd', 'compressed')
                self._rhofile = File(s + '_rhoOpt.pvd', 'compressed')
                self._pfile = File(s + '_pOpt.pvd', 'compressed')
            else:
                self._ufile = File(s + '_u.pvd', 'compressed')
                self._rhofile = File(s + '_rho.pvd', 'compressed')
                self._pfile = File(s + '_p.pvd', 'compressed')
            self._uDualfile = File(s + '_uDual.pvd', 'compressed')
            self._rhoDualfile = File(s + '_rhoDual.pvd', 'compressed')
            self._pDualfile = File(s + '_pDual.pvd', 'compressed')
            self.meshfile = File(s + '_mesh.xml')
        else:
            self._ufile = File(s + '_u%d.pvd' % n, 'compressed')
            self._rhofile = File(s + '_rho%d.pvd' % n, 'compressed')
            self._pfile = File(s + '_p%d.pvd' % n, 'compressed')
            self._uDualfile = File(s + '_uDual%d.pvd' % n, 'compressed')
            self._rhoDualfile = File(s + '_rhoDual%d.pvd' % n, 'compressed')
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
                self._rhofile << tau
            else:
                self._uDualfile << u
                self._pDualfile << p
                self._rhoDualfile << tau

    def Plot(self, problem, W, w):
        u = w.split()[0]
        tau = w.split()[2]
        p = w.split()[1]

        if self.vizU is None:
            # Plot velocity and pressure
            self.vizU = plot(u, title='Velocity', rescale=True)
            self.vizP = plot(p, title='Pressure', rescale=True, elevate=0.0)
            self.vizRho = plot(tau, title='Stress', rescale=True, elevate=0.0)
        else:
            self.vizU.plot(u)
            self.vizP.plot(p)
            self.vizRho.plot(tau)

    def __str__(self):
        return 'OldroydB'
