__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"

from solverbase import *

class Solver(SolverBase):
    '''
        Solver class for Peregrine system.
    '''

    def __init__(self, options):
        SolverBase.__init__(self, options)

        #Parameters
        self.Re = None
        self.Ro = None
        self.Hfile = None
        self.vizH = None

    #define functions spaces
    def function_space(self, mesh):
        '''
            The Peregrine system doesn't quite fit into our normal scheme, so we
            need to define the functions space for the Peregrine System. This is
            mostly so that the wave object is represented.
        '''
        V = VectorFunctionSpace(mesh, 'CG', self.Pu)
        Q = FunctionSpace(mesh, 'CG', self.Pp)
        self.Q = Q #Bad hack for being able to project into Q space
        W = MixedFunctionSpace([V, Q])

        self.Zeta = Function(Q, name='zeta')
        self.Zeta_ = Function(Q, name='zeta_')
        self.Zeta__ = Function(Q, name='zeta__')
        self.H = Function(Q, name='Bathymetry')
        self.H_ = Function(Q, name='PreviousBathymetry')

        return W

    def strong_residual(self, problem, H, U, v, eta, chi):
        '''
            Defines the strong residual for Peregrine System.
        '''
        #Parameters
        sigma = problem.sigma
        epsilon = problem.epsilon

        k = problem.k

        zeta_tt = 1./k**2*(self.Zeta - 2*self.Zeta_ + self.Zeta__)
        zeta_t = 1./k*(self.Zeta - self.Zeta_)

        #strong form for stabilization
        R1 = epsilon*grad(v)*U + grad(chi)
        z1 = -sigma**2*H/2.*grad(zeta_tt)

        R2 = div(v*(H + epsilon*eta))
        z2 = zeta_t

        return R1, R2, z1, z2

    def weak_residual(self, problem, k, W, w, w_, wt, ei_mode=False):
        '''
            Defines the weak residual for Peregrine System, including LS
            Stabilization.
        '''
        (U, eta) = (as_vector((w[0], w[1])), w[2])
        (U_, eta_) = (as_vector((w_[0], w_[1])), w_[2])
        (v, chi) = (as_vector((wt[0], wt[1])), wt[2])

        h = CellSize(W.mesh()) #mesh size
        d = 0.1*h**(3./2.) #stabilization parameter
        d1, d2 = self.stabilization_parameters(U_, eta_, k, h) #stabilization parameters

        #set up error indicators
        Z = FunctionSpace(W.mesh(), "DG", 0)
        z = TestFunction(Z)

        alpha = self.alpha #time stepping method

        #Parameters
        sigma = problem.sigma
        epsilon = problem.epsilon

        t0 = problem.t0 #initial time
        self.t0 = t0
        T = problem.T #Final time

        #set up our wave object
        self.wave_object(problem, self.Q, t0, k)

        zeta_tt = 1./k**2*(self.Zeta - 2*self.Zeta_ + self.Zeta__)
        zeta_t = 1./k*(self.Zeta - self.Zeta_)

        #Time stepping method
        U_alpha = (1. - alpha)*U_ + alpha*U
        eta_alpha = (1. - alpha)*eta_ + alpha*eta
        zeta_alpha = (1. - alpha)*self.Zeta_ + alpha*self.Zeta
        H_alpha = (1. - alpha)*self.H_ + alpha*self.H

        t = t0 + k

        if(not self.options["stabilize"] or ei_mode):
          d = 0
          d1 = 0
          d2 = 0
        if(not ei_mode):
          z = 1.

        #weak form of the equations
        r = z*(1./k*inner(U-U_,v) + epsilon*inner(grad(U_alpha)*U_alpha,v) \
            - div(v)*eta_alpha)*dx

        r += z*(sigma**2*1./k*div(H_alpha*(U-U_))*div(H_alpha*v/2.) \
              - sigma**2*1./k*div(U-U_)*div(H_alpha**2*v/6.))*dx
        r += z*sigma**2*zeta_tt*div(H_alpha*v/2.)*dx

        r += z*(1./k*(eta-eta_)*chi + zeta_t*chi)*dx
        r -= z*inner(U_alpha,grad(chi))*(epsilon*eta_alpha + H_alpha)*dx

        r += z*d*(inner(grad(U_alpha),grad(v)) + inner(grad(eta_alpha),grad(chi)))*dx

        R1, R2, z1, z2 = self.strong_residual(problem, H_alpha, U_alpha, U_alpha, \
                eta_alpha, eta_alpha)
        Rv1, Rv2, z1, z2 = self.strong_residual(problem, H_alpha, U_alpha, v, \
                eta_alpha, chi)
        r += z*(d1*inner(R1 + z1, Rv1) + d2*(R2 + z2)*Rv2)*dx

        return r

    def stabilization_parameters(self, U_, eta_, k, h):
        K1  = 2.
        K2  = 2.
        d1 = K1*(k**(-2) + inner(U_,U_)*h**(-2))**(-0.5)
        d2 = K2*(k**(-2) + eta_*eta_*h**(-2))**(-0.5)

        return d1, d2

    #Optimization Function
    def Optimize(self, problem, w):
        '''
            Shape optimization for Peregrine System.
        '''
        eta = w.split()[1]
        self.optfile = File(self.s + '_Opt.pvd') #file for solution to optimization

        #bounds on object
#        lb = project(Expression('-0.5'), self.Q, name='LowerBound')
#        ub = project(Expression('0.0'), self.Q, name='UpperBound')

        #Functionnal to be minimized: L2 norm over a subdomain
        J = Functional(-inner(eta, eta)*dx*dt[FINISH_TIME] + self.Zeta*self.Zeta*dx)

        #shape parameters
        m = [Control(p) for p in problem.params]
        Jhat = ReducedFunctional(J, m) #Reduced Functional
        opt_params = minimize(Jhat)
        opt_object = project(problem.object_init(opt_params), self.Q)
        if self.options['plot_solution']:
            plot(opt_object, title='Optimization result.')
            interactive()
        else:
            optfile << opt_object

        print 'H=%f, N1=%f, N2=%f, W=[%f, %f, %f, %f], dz=%f]' % \
                (opt_params[0], opt_params[1], opt_params[2], opt_params[3], \
                opt_params[4], opt_params[5], opt_params[6], opt_params[7])

    def functional(self,mesh,w):
        '''
            Functional for mesh adaptivity
        '''

        (u, eta) = (as_vector((w[0], w[1])), w[2])

        M = u[0]*dx # Mean of the x-velocity in the whole domain

        return M

    def wave_object(self, problem, Q, t, k):
        '''
            Update for wave object at each time step.
        '''
        H, H_, zeta, zeta_, zeta__ = problem.update_bathymetry(Q, t)

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

        if self.options['save_frequency'] !=0 and (self._timestep - 1) % self.options['save_frequency'] == 0:
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
        if n==-1:
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
        else :
            self.vizU.plot(u)
            self.vizP.plot(eta)
            self.vizH.plot(self.H)

    def __str__(self):
          return 'Peregrine'
