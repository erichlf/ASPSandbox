from AFES import *
import sys
from time import time

# Constants related to the geometry
bmarg = 1.e-3 + DOLFIN_EPS
xmin = 0.0
xmax = 2.2
ymin = 0.0
ymax = 0.41
xcenter = 0.2
ycenter = 0.2
radius = 0.05
Diameter = 2. * radius
Um = 1.5  # max velocity

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


class NoSlipBoundary(SubDomain):

    def __init__(self):
        SubDomain.__init__(self)
        self.object = ObjectBoundary()

    def inside(self, x, on_boundary):

        return on_boundary and (near(x[1], ymin) or near(x[1], ymax)
                                or self.object.inside(x, on_boundary))


class ObjectBoundary(SubDomain):

    def __init__(self):
        SubDomain.__init__(self)

    def inside(self, x, on_boundary):
        dx = x[0] - xcenter
        dy = x[1] - ycenter
        r = sqrt(dx * dx + dy * dy)

        return on_boundary and r < radius + bmarg


class OutflowBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], xmax)


def initial_conditions(self, W):  # create initial conditions

    return project(InitialConditions(), W)


def time_step(self, u, mesh):
    C_CFL = 100.
    return C_CFL * mesh.hmin() / u


def boundary_conditions(self, W, t):
    self.U.t = t
    # Create inflow boundary condition
    bc0 = DirichletBC(W.sub(0), self.U, InflowBoundary())

    # Create no-slip boundary condition
    bc1 = DirichletBC(W.sub(0), Constant((0, 0)), NoSlipBoundary())

    # Create outflow boundary condition for pressure
    bc2 = DirichletBC(W.sub(1), Constant(0), OutflowBoundary())

    # Collect boundary conditions
    bcs = [bc0, bc1, bc2]

    return bcs


def F(self, t):  # forcing function for the momentum equation
    return Expression(('0', '0'))


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


def strong_residual(W, w, w2):  # strong residual for cG(1)cG(1)
    (u, p) = (as_vector((w[0], w[1])), w[2])
    (U, P) = (as_vector((w2[0], w2[1])), w2[2])

    R1 = grad(U) * u + grad(p)
    R2 = div(u)

    return R1, R2


def weak_residual(self, problem, k, W, w, ww, w_, wt, ei_mode=False):

    (u, p) = (as_vector((w[0], w[1])), w[2])
    (U, P) = (as_vector((ww[0], ww[1])), ww[2])
    (U_, P_) = (as_vector((w_[0], w_[1])), w_[2])
    (v, q) = (as_vector((wt[0], wt[1])), wt[2])

    h = CellSize(W.mesh())  # mesh size

    nu = problem.nu  # Reynolds Number

    t0 = problem.t0

    t = t0 + k
    # forcing and mass source/sink
    F = problem.F(t)

    d1, d2 = stabilization_parameters(h)

    # weak residual for cG(1)cG(1)
    r = (1. / k) * inner(U - U_, v) * dx \
        + inner(grad(p) + grad(u) * u, v) * dx
    r += nu * inner(grad(u), grad(v)) * dx
    r += div(u) * q * dx

    # forcing function
    r -= inner(F, v) * dx

    # least squares stabilization
    R1, R2 = strong_residual(W, w, w)
    Rv1, Rv2 = strong_residual(W, wt, w)
    r += (d1 * inner(R1 - F, Rv1) + d2 * R2 * Rv2) * dx

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

    problem.nu = Constant(1E-3)

    # Spacial discretization
    Nx = 60

    # setup our domain and conditions
    channel = Rectangle(xmin, ymin, xmax, ymax)
    bluff = Circle(xcenter, ycenter, radius)
    domain = channel - bluff
    problem.mesh = Mesh(domain, Nx)

    H = ymax
    problem.U = Expression(('4*Um*x[1]*(H - x[1])/(H*H)*t*t/(1+t*t)',
                            '0.0'), Um=Um, H=ymax, t=problem.t0)
    problem.Ubar = 4. / 3. * Um * ymax * (H - ymax / 2.) / (H * H)

    # Reynolds number
    problem.Re = problem.Ubar * Diameter * problem.nu

    # set timestep based on UFL type condition
    Problem.time_step = time_step
    problem.k = problem.time_step(problem.Ubar, problem.mesh)

    # create initial and boundary conditions
    Problem.initial_conditions = initial_conditions
    Problem.boundary_conditions = boundary_conditions

    # create some problem specific functions
    Problem.F = F
    Problem.update = update_problem

    # create solver instance
    solver = Solver(options)

    # create our user defined problem and function space
    Solver.function_space = function_space
    Solver.weak_residual = weak_residual

    Solver.__str__ = 'NSE'
    Problem.__str__ = 'Cylinder'

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
