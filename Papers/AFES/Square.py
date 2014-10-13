from problembase import *


class Problem(ProblemBase):

    def __init__(self, options):  # initialize the problem
        ProblemBase.__init__(self, options)  # send options to problembase

        Nx = options['Nx']
        Ny = options['Ny']

        # Create mesh
        self.mesh = UnitSquareMesh(20, 20)

        # time domain and time step
        self.t0 = 0.
        self.T = 5.
        self.k = 0.01

        self.kappa = Constant(1E-2)  # heat capacity

        self.u0 = Expression('sin(pi*x[0])*sin(pi*x[1])')  # ic
        self.f = Constant(0)  # forcing
        self.g = Constant(0)  # boundary condition

    # create initial conditions
    def initial_conditions(self, V):

        return project(self.u0, V)

    # Create boundary condition
    def boundary_conditions(self, V, t):

        return DirichletBC(V, self.g, 'on_boundary')

    def __str__(self):
        return 'Square'
