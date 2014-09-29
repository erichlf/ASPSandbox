__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"

from solverbase import *

class Solver(SolverBase):
    '''
        Solver class for solving the variable density incompressible
        Navier-Stokes equation.
    '''

    def __init__(self, options):
        SolverBase.__init__(self, options)

        self.Re = options['Re']
        self.Fr = options['Fr']
        self.vizRho = None
        self.rhofile = None

    # Define function spaces
    def function_space(self, mesh):
        V = VectorFunctionSpace(mesh, 'CG', self.Pu)
        R = FunctionSpace(mesh, 'CG', self.Pp)
        Q = FunctionSpace(mesh, 'CG', self.Pp)

        W = MixedFunctionSpace([V, R, Q])

        return W

    #strong residual for cG(1)cG(1)
    def strong_residual(self, w, w2):# U, v, rho, r, p):
        (u, rho, p) = (as_vector((w[0], w[1])), w[2], w[3])
        (U, Rho, P) = (as_vector((w2[0], w2[1])), w2[2], w2[3])

        R1 = Rho*grad(U)*u + grad(p)
        R2 = div(Rho*u)
        R3 = div(u)

        return R1, R2, R3

    #weak residual for cG(1)cG(1)
    def weak_residual(self, problem, k, W, w, ww, w_, wt, ei_mode=False):
        #rho = 1/rho' - 1
        (u, rho, p) = (as_vector((w[0], w[1])), w[2], w[3])
        (U, Rho, P) = (as_vector((ww[0], ww[1])), ww[2], ww[3])
        (U_, Rho_, P_) = (as_vector((w_[0], w_[1])), w_[2], w_[3])
        (v, nu, q) = (as_vector((wt[0], wt[1])), wt[2], wt[3])

        h = CellSize(W.mesh()) #mesh size
        d1, d2, d3 = self.stabilization_parameters(U_, Rho_, P_, k, h) #stabilization parameters
        d = 0.1*h**(3./2.) #stabilization parameter

        #set up error indicators
        Z = FunctionSpace(W.mesh(), "DG", 0)
        z = TestFunction(Z)

        Re = self.Re #Reynolds Number

        t0 = problem.t0

        t = t0 + k
        #forcing and mass source/sink
        f = problem.F1(t)

        #least squares stabilization
        if(not self.options["stabilize"] or ei_mode):
          d1 = 0
          d2 = 0
          d3 = 0
          d = 0
        if(not ei_mode):
          z = 1.

        #weak form of the equations
        r = z*((1./k)*(Rho - Rho_) + div(rho*u))*nu*dx #mass equation
        r += z*(rho*((1./k)*inner(U - U_,v) \
            + inner(grad(u)*u,v)) \
            + inner(grad(p),v) \
            + 1./Re*inner(grad(u),grad(v)) \
            + rho*dot(f,v))*dx #momentum equation
        r += z*div(u)*q*dx #continuity equation

        R1, R2, R3 = self.strong_residual(w,w)
        Rv1, Rv2, Rv3 = self.strong_residual(wt,w)
        r += z*(d1*inner(R1 - rho*f, Rv1) + d2*R2*Rv2 + d3*R3*Rv3)*dx
        r += z*d*(inner(grad(u),grad(v)))*dx
        r += z*d*(inner(grad(rho),grad(nu)))*dx

        return r

    def stabilization_parameters(self, u, rho, p, k, h):
        K1  = 1.
        K2  = 1.
        K3  = 0.5
        d1 = K1*(k**(-2) + inner(u,u)*h**(-2))**(-0.5)
        d2 = K2*(k**(-2) + rho*rho*h**(-2))**(-0.5)
        d3 = K3*h

        return d1, d2, d3

    def Save(self, problem, w, dual=False):
        u = w.split()[0]
        rho = w.split()[1]
        p = w.split()[2]

        if self.options['save_frequency'] !=0 and (self._timestep - 1) % self.options['save_frequency'] == 0:
            if not dual:
                self._ufile << u
                self._pfile << p
                self._rhofile << rho
            else:
                self._uDualfile << u
                self._pDualfile << p
                self._rhoDualfile << rho

    def file_naming(self, n=-1, dual=False):
        if n==-1:
            self._ufile = File(self.s + '_u.pvd', 'compressed')
            self._rhofile = File(self.s + '_rho.pvd', 'compressed')
            self._pfile = File(self.s + '_p.pvd', 'compressed')
            self._uDualfile = File(self.s + '_uDual.pvd', 'compressed')
            self._rhoDualfile = File(self.s + '_rhoDual.pvd', 'compressed')
            self._pDualfile = File(self.s + '_pDual.pvd', 'compressed')
        else:
            self._ufile = File(self.s + '_u%d.pvd' % n, 'compressed')
            self._rhofile = File(self.s + '_rho%d.pvd' % n, 'compressed')
            self._pfile = File(self.s + '_p%d.pvd' % n, 'compressed')
            self._uDualfile = File(self.s + '_uDual%d.pvd' % n, 'compressed')
            self._rhoDualfile = File(self.s + '_rhoDual%d.pvd' % n, 'compressed')
            self._pDualfile = File(self.s + '_pDual%d.pvd' % n, 'compressed')

    def Plot(self, problem, W, w):
        u = w.split()[0]
        rho = w.split()[1]
        p = w.split()[2]

        if self.vizU is None:
            # Plot velocity and pressure
            self.vizU = plot(u, title='Velocity', rescale=True)
            self.vizP = plot(p, title='Pressure', rescale=True, elevate=0.0)
            self.vizRho = plot(rho, title='Density', rescale=True, elevate=0.0)
        else :
            self.vizU.plot(u)
            self.vizP.plot(p)
            self.vizRho.plot(rho)

    def __str__(self):
          return 'DensityNSE'
