from AFES import *
import sys
import time


'''
    This class is used to define our domain, bcs, and ics.
'''


def initial_conditions(self, V):  # create initial conditions

    return project(self.u0, V)


def boundary_conditions(self, V, t):  # Create boundary condition

    return DirichletBC(V, self.g, 'on_boundary')


'''
    The following will create our solver and tells AFES what to do with
    the solutions that you obtain.
'''


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


def main():
    # setup problem options
    options = {
        'T': 1.,  # stopping time
        'k': 0.01,  # time step
        'theta': 1,  # implicit Euler
        'kappa': Constant(1E-2),  # heat coefficient
        'u0': Expression('sin(pi*x[0])*sin(pi*x[1])'),  # ic
        'f': Expression('0'),  # forcing
        'g': Expression('0'),  # boundary condition
        'Nx': 20,  # spacial discretization along x-axis
        'Ny': 20,  # spacial discretization along y-axis
        'velocity_order': 1,  # order of the finite element

        # setup AFES specific options
        'adaptive': False,  # mesh adaptivity
        'optimize': False,  # optimize as defined in solver
        'save_solution': False,  # don't save the solution
        'plot_solution': True,  # plot the solution
    }

    problem = Problem(options)

    # Spacial discretization
    Nx = options['Nx']
    Ny = options['Ny']

    problem.mesh = UnitSquareMesh(Nx, Ny)  # Create mesh

    # time domain and time step
    problem.t0 = 0.
    problem.T = options['T']
    problem.k = options['k']

    problem.kappa = options['kappa']  # heat capacity

    problem.u0 = options['u0']  # ic
    problem.f = options['f']  # forcing
    problem.g = options['g']  # boundary condition

    Problem.initial_conditions = initial_conditions  # create initial conditions
    Problem.boundary_conditions = boundary_conditions  # Create boundary condition

    # Create problem and solver
    solver = Solver(options)

    # create our user defined problem and function space
    Solver.function_space = function_space
    Solver.weak_residual = weak_residual
    Solver.Plot = Plot  # Add the plotter

    Solver.__str__ = 'Heat'
    Problem.__str__ = 'Square'

    # Solve problem with solver
    wct = time.time()
    w = solver.solve(problem)

    # Compute elapsed time
    wct = time.time() - wct

    sys.stdout.flush()
    sys.stdout.write('\033[K')
    print 'Solved %s for the %s problem in %g seconds.' % (solver.__str__,
                                                           problem.__str__, wct)

    return 0

if __name__ == "__main__":
    sys.exit(main())
