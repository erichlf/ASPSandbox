from AFES import *
import sys
from time import time

# Constants related to the geometry
xmin = 0.
xmax = 17.
ymin = 0.
ymax = 2.
Ud = Expression(('0', '0'))  # domain velocity for ALE terms
n = 3  # nth resonance
c = 350.  # speed of sound in m/s
rho = 1.225  # density of air in kg/m^3
f = (2 * n - 1) * c / (4 * (xmax - xmin))  # nth resonance frequency

'''
    The following functions are used to define ics, and bcs.
'''


class InitialConditions(Expression):

    def eval(self, values, x):
        values[0] = 0.
        values[1] = 0.
        values[2] = 0.

    def value_shape(self):
        return (3,)


class InflowBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], xmin)


class NoNormalBoundary(SubDomain):

    def __init__(self):
        SubDomain.__init__(self)

    def inside(self, x, on_boundary):

        return on_boundary and (near(x[1], ymin) or near(x[1], ymax))


class OutflowBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], xmax)


def initial_conditions(self, W):  # create initial conditions

    return project(InitialConditions(), W)


def boundary_conditions(self, W, t):
    self.Ug.t = t
    # Create inflow boundary condition
    bc0 = DirichletBC(W.sub(0), self.Ug, InflowBoundary())

    # Create no-normal flow boundary condition
    bc1 = DirichletBC(W.sub(0).sub(1), Constant(0), NoNormalBoundary())

    # Create outflow boundary condition for pressure
    bc2 = DirichletBC(W.sub(1), Constant(0), OutflowBoundary())

    # Collect boundary conditions
    bcs = [bc0, bc1, bc2]

    return bcs


def F(self, t):  # forcing function for the momentum equation
    return Expression(('0', '0'))


def Q(self, t):  # forcing function for the continuity equation
    return Expression('0')


def update_problem(self, W, t):

    bcs = self.boundary_conditions(W, t)

    return bcs


'''
    The following will create our solver and tells AFES what to do with
    the solutions that you obtain.
'''


def function_space(self, mesh):
    V = VectorFunctionSpace(mesh, 'CG', 1)
    Q = FunctionSpace(mesh, 'CG', 1)
    W = MixedFunctionSpace([V, Q])

    return W


def strong_residual(self, W, w, w2):  # strong residual for cG(1)cG(1)
    (u, p) = (as_vector((w[0], w[1])), w[2])

    rho = self.rho  # density
    c = self.c

    R1 = -rho * grad(u) * Ud + grad(p)
    R2 = -1. / (rho * c * c) * dot(Ud, grad(p)) + div(u)

    return R1, R2


def weak_residual(self, problem, k, W, w, ww, w_, wt, ei_mode=False):

    (u, p) = (as_vector((w[0], w[1])), w[2])
    (U, P) = (as_vector((ww[0], ww[1])), ww[2])
    (U_, P_) = (as_vector((w_[0], w_[1])), w_[2])
    (v, q) = (as_vector((wt[0], wt[1])), wt[2])

    h = CellSize(W.mesh())  # mesh size

    rho = problem.rho  # density
    c = problem.c  # speed of sound
    self.rho = rho
    self.c = c

    t0 = problem.t0

    t = t0 + k
    # forcing and mass source/sink
    F = problem.F(t)
    Q = problem.Q(t)

    d1, d2 = stabilization_parameters(h)

    # weak residual for cG(1)cG(1)
    r = (rho / k) * inner(U - U_, v) * dx \
        + inner(grad(p) + rho * grad(u) * Ud, v) * dx
    r += 1. / (rho * c * c) * (P - P_) * q * dx \
        - (1. / (rho * c * c) * dot(Ud, grad(p)) - div(u)) * q * dx

    # forcing function
    r -= (inner(F, v) + Q * q) * dx

    # least squares stabilization
    R1, R2 = self.strong_residual(W, w, w)
    Rv1, Rv2 = self.strong_residual(W, wt, w)
    r += (d1 * inner(R1 - F, Rv1) + d2 * (R2 - Q) * Rv2) * dx

    return r


def stabilization_parameters(h):
    K1 = 1.
    K2 = 1.
    d1 = K1 * h
    d2 = K2 * h

    return d1, d2


def main():
    # setup problem options
    options = OPTIONS.copy()
    options['T'] = 10.

    # create problem instance
    problem = Problem(options)

    # Spacial discretization
    Nx = 10
    Ny = 100

    # setup our domain
    problem.mesh = RectangleMesh(xmin, ymin, xmax, ymax, Nx, Ny)

    problem.Ug = Expression(('sin(2*pi*f*t)', '0.0'), f=f, t=problem.t0)

    # problem specific parameters
    problem.rho = rho
    problem.c = c

    problem.k = 0.01  # time step

    # create initial and boundary conditions
    Problem.initial_conditions = initial_conditions
    Problem.boundary_conditions = boundary_conditions

    # create some problem specific functions
    Problem.F = F
    Problem.Q = Q
    Problem.update = update_problem

    # create solver instance
    solver = Solver(options)

    # create our user defined problem and function space
    Solver.function_space = function_space
    Solver.weak_residual = weak_residual
    Solver.strong_residual = strong_residual

    Solver.__str__ = 'MixedWave'
    Problem.__str__ = 'Channel'

    # Solve problem with solver
    wct = time()
    w = solver.solve(problem)
    wct = time() - wct  # Compute elapsed time

    sys.stdout.flush()
    sys.stdout.write('\033[K')
    print 'Solved %s for the %s problem in %g seconds.' % (solver.__str__,
                                                           problem.__str__, wct)

    return 0

if __name__ == "__main__":
    sys.exit(main())
