__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2014-08-12"
__copyright__ = "Copyright (C) 2013-2014 " + __author__
__license__ = "GNU GPL version 3 or any later version"

from ASP import *
from ASP import Problem as ProblemBase
from mshr import *

# Constants related to the geometry
bmarg = 1.E-3 + DOLFIN_EPS
xmin, xmax, ymin, ymax = 0.0, 2.2, 0.0, 0.41
xcenter, ycenter = 0.2, 0.2
radius = 0.05
Diameter = 2. * radius
Um = 1.5  # max velocity


class InitialConditions(Expression):

    def eval(self, values, x):
        values[0:3] = 0., 0., 0.

    def value_shape(self):
        return (3,)


class InflowBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], xmin)


class NoSlipBoundary(SubDomain):

    def inside(self, x, on_boundary):
        dx, dy = x[0] - xcenter, x[1] - ycenter
        r = sqrt(dx * dx + dy * dy)

        return on_boundary and (near(x[1], ymin) or near(x[1], ymax)
                                or r < radius + bmarg)


class OutflowBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], xmax)


class Problem(ProblemBase):

    '''
        Provides the 2D flow around a cylinder body. This problem is described
        in "Benchmark Computations of Laminar Flow Around a Cylinder" by M.
        Schaefer and S. Turek.
    '''

    def __init__(self, options, cube=False):

        global xmax, xcenter, ycenter

        # Load mesh
        options['Ny'], options['Nz'] = None, None
        ProblemBase.__init__(self, options)
        self.Nx = options['Nx']

        try:
            self.nu = options['nu']
        except:
            self.nu = 1E-3

        self.t0 = 0.
        try:
            self.T = options['T']
        except:
            self.T = 10

        H = ymax
        # setup our domain and conditions
        try:
            self.mesh = Mesh(options['initial_mesh'])
        except:
            channel = Rectangle(Point(xmin, ymin), Point(xmax, ymax))
            bluff = Circle(Point(xcenter, ycenter), radius)
            domain = channel - bluff
            self.mesh = generate_mesh(domain, self.Nx)

        self.noSlip = Constant((0, 0))
        self.U = Expression(('4*Um*x[1]*(H - x[1])/(H*H)*t*t/(1+t*t)',
                            '0.0'), Um=Um, H=ymax, t=self.t0)
        self.Ubar = 4. / 3. * Um * ymax * (H - ymax / 2.) / (H * H)

        self.k = self.time_step(self.Ubar, self.mesh)  # mesh size

    def initial_conditions(self, W):
        w0 = InitialConditions()

        w0 = project(w0, W)

        return w0

    def time_step(self, u, mesh):
        C_CFL = 100.
        return C_CFL * mesh.hmin()/u

    def boundary_conditions(self, W, t):
        self.U.t = t
        # Create inflow boundary condition
        bc0 = DirichletBC(W.sub(0), self.U, InflowBoundary())

        # Create no-slip boundary condition
        bc1 = DirichletBC(W.sub(0), self.noSlip,
                          NoSlipBoundary(self.dim, self.cube))

        # Create outflow boundary condition for pressure
        bc2 = DirichletBC(W.sub(1), Constant(0), OutflowBoundary())

        bcs = [bc0, bc1, bc2]  # Collect boundary conditions

        return bcs

    def F1(self, t):  # forcing function for the momentum equation
        f = Constant((0, 0))

        return f

    def F2(self, t):  # mass source for the continuity equation
        return Constant(0)

    def update(self, W, t):  # update the bc for each time step

        return self.boundary_conditions(W, t)

    def functional(self, W, w):
        '''
            This is the functional used for adaptivity.
            We assume the problem is much like NSE.
        '''
        u = as_vector((w[0], w[1]))

        return u[0] * dx  # Mean of the x-velocity in the whole domain

    def __str__(self):
        return 'Cylinder'
