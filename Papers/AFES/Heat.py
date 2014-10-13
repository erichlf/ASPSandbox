from solverbase import *


class Solver(SolverBase):

    def __init__(self, options):  # initialize the solver
        options['theta'] = 1.  # Choose implicit Euler
        SolverBase.__init__(self, options)  # send options to solverbase

    def function_space(self, mesh):
        V = FunctionSpace(mesh, 'CG', 1)  # Define function spaces

        return V

    # define our weak form
    def weak_residual(self, problem, k, V, u, U, U_, v, ei_mode=False):
        kappa = problem.kappa

        f = problem.f  # forcing and mass source/sink

        # weak form of the equations
        r = (1. / k) * (U - U_) * v * dx \
            + kappa * inner(grad(u), grad(v)) * dx \
            - f * v * dx

        return r

    def __str__(self):
        return 'Heat'

    def Plot(self, problem, V, u):

        # Plot energy
        if self.vizU is None:
            self.vizU = plot(u, title='Temperature', rescale=False)
        else:
            self.vizU.plot(u)

    def Save(self, problem, u, dual=False):

        self._ufile << u

    def file_naming(self, n=-1, dual=False):
        s = 'HeatEq'

        # name our files
        self._ufile = File(s + '_u.pvd', 'compressed')
        self.meshfile = File(s + '_mesh.xml')
