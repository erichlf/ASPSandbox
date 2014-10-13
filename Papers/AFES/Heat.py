from AFES import *


class Problem(ProblemBase):

    '''
        This class is used to define our domain, bcs, and ics.
    '''

    def __init__(self, options):  # initialize the problem
        ProblemBase.__init__(self, options)  # send options to problembase

        # Spacial discretization
        Nx = options['Nx']
        Ny = options['Ny']

        self.mesh = UnitSquareMesh(Nx, Ny)  # Create mesh

        # time domain and time step
        self.t0 = 0.
        self.T = options['T']
        self.k = options['k']

        self.kappa = options['kappa']  # heat capacity

        self.u0 = options['u0']  # ic
        self.f = options['f']  # forcing
        self.g = options['g']  # boundary condition

    # create initial conditions
    def initial_conditions(self, V):

        return project(self.u0, V)

    # Create boundary condition
    def boundary_conditions(self, V, t):

        return DirichletBC(V, self.g, 'on_boundary')

    def __str__(self):
        return 'Square'


class Solver(SolverBase):

    '''
        This class is used to create our solver and tells AFES what to do with
        the solutions that you obtain.
    '''

    def __init__(self, options):  # initialize the solver
        SolverBase.__init__(self, options)  # send options to solverbase

    def function_space(self, mesh):
        V = FunctionSpace(mesh, 'CG', 1)  # Define function spaces

        return V

    # define our weak form
    def weak_residual(self, problem, k, V, u, U, U_, v, ei_mode=False):
        kappa = problem.kappa  # heat coefficient

        f = problem.f  # source/sink

        # weak form of the equations
        r = (1. / k) * (U - U_) * v * dx \
            + kappa * inner(grad(u), grad(v)) * dx \
            - f * v * dx

        return r

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

    def __str__(self):
        return 'Heat'


def main():
    # setup problem options
    options['T'] = 1.  # stopping time
    options['k'] = 0.01  # time step
    options['theta'] = 1  # implicit Euler
    options['kappa'] = Constant(1E-2)  # heat coefficient
    options['u0'] = Expression('sin(pi*x[0])*sin(pi*x[1])')  # ic
    options['f'] = Expression('0')  # forcing
    options['g'] = Expression('0')  # boundary condition

    # setup AFES specific options
    options['adaptive'] = False  # mesh adaptivity
    options['optimize'] = False  # optimize as defined in solver
    options['save_solution'] = False  # don't save the solution
    options['plot_solution'] = True  # plot the solution

    # Create problem and solver
    solver = Solver(options)
    problem = Problem(options)

    solver_name = solver.__str__()
    problem_name = problem.__str__()

    # Solve problem with solver
    wct = time.time()
    w = solver.solve(problem)

    # Compute elapsed time
    wct = time.time() - wct

    sys.stdout.flush()
    sys.stdout.write('\033[K')
    print 'Solved %s for the %s problem in %g seconds.' % (solver_name,
                                                           problem_name, wct)

    return 0

if __name__ == "__main__":
    sys.exit(main())
