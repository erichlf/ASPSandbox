__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__  = "GNU GPL version 3 or any later version"

from solverbase import *

class Solver(SolverBase):
    '''
        Solver class for solving the Advection-Diffusion-Reaction equation.
    '''

    def __init__(self, options):
        SolverBase.__init__(self, options)
        self.epsilon = options['epsilon']
        self.a = options['a']
        self.beta = options['beta']

    #strong residual for cG(1)cG(1)
    def strong_residual(self,u):
        epsilon = self.epsilon #diffusion coefficient
        a = self.a #reaction coefficient
        beta = self.beta #velocity

        R = a*u + dot(beta,grad(u))

        return R

    #weak residual for cG(1)cG(1)
    def weak_residual(self, problem, W, w, w_, wt, ei_mode=False):
        h = CellSize(W.mesh()) #mesh size
        d = self.stabilization_parameters(U_,h) #stabilization parameters

        #set up error indicators
        Z = FunctionSpace(W.mesh(), "DG", 0)
        z = TestFunction(Z)

        a = self.a #reaction coefficient
        epsilon = self.epsilon #diffusion coefficient
        beta = self.beta #velocity

        alpha = self.alpha #time stepping method
        self.k = Expression('dt', dt=problem.k)
        t0 = problem.t0

        #U_(k+alpha)
        U_alpha = (1.0-alpha)*U_ + alpha*U

        t = t0 + self.k
        #forcing and mass source/sink
        F = problem.F1(t)

        #least squares stabilization
        if(not self.options["stabilize"] or ei_mode):
          d = 0
        if(not ei_mode):
          z = 1.

        #weak form of the equations
        r = z*((1./self.k)*inner(U - U_,v) \
            + alpha*inner(u,v)
            + epsilon*inner(grad(U_alpha),grad(v)) \
            + inner(dot(beta,grad(U_alpha)),v)*dx

        #forcing function
        r -= z*inner(F,v)*dx

        R = self.strong_residual(U_alpha)
        Rv = self.strong_residual(v)
        r += z*d1*inner(R1 - F, Rv1)*dx

        return r

    def stabilization_parameters(self,U,h):
        K  = 1.
        if(h > self.epsilon):
            d = K*(self.k**(-2) + inner(U,U)*h**(-4))**(-0.5)
        else:
            d = K*(self.k**(-2) + inner(U,U)*h**(-2))**(-0.5)

        return d

    #this is the functional used for adaptivity
    def functional(self, mesh, w):
        M = u[0]*dx # Mean of the x-velocity in the whole domain

        return M

    def __str__(self):
          return 'ADR'
