__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2014-08-12"
__copyright__ = "Copyright (C) 2013-2014 " + __author__
__license__ = "GNU GPL version 3 or any later version"
#
#   adapted from ns in nsbench originally developed by
#   Anders Logg <logg@simula.no>
#

from AFES import *
from AFES import Problem as ProblemBase
from mshr import *

# Constants related to the geometry
bmarg = 1.e-3 + DOLFIN_EPS
xmin = 0.0
xmax = 2.2
ymin = 0.0
ymax = 0.41
zmin = 0.0
zmax = 0.41
xcenter = 0.2
ycenter = 0.2
radius = 0.05
Diameter = 2. * radius
Um = 1.5  # max velocity


class InitialConditions2D(Expression):

    def eval(self, values, x):
        values[0] = 0.
        values[1] = 0.
        values[2] = 0.
        values[3] = 0.
        values[4] = 0.
        values[5] = 0.
        values[6] = 0.

    def value_shape(self):
        return (7,)


class InitialConditions3D(Expression):

    def eval(self, values, x):
        values[0] = 0.
        values[1] = 0.
        values[2] = 0.
        values[3] = 0.

    def value_shape(self):
        return (4,)


# Inflow boundary
class InflowBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], xmin)


# No-slip boundary
class NoSlipBoundary(SubDomain):

    def __init__(self, dim,):
        SubDomain.__init__(self)
        self.dim = dim
        self.object = ObjectBoundary(dim)

    def inside(self, x, on_boundary):

        return on_boundary and (near(x[1], ymin) or near(x[1], ymax)
                                or (self.dim == 3 and (near(x[2], zmin)
                                                       or near(x[2], zmax)))
                                or self.object.inside(x, on_boundary))


class ObjectBoundary(SubDomain):

    def __init__(self, dim):
        SubDomain.__init__(self)
        self.dim = dim

    def inside(self, x, on_boundary):

        return on_boundary and (x[0] > xmin + DOLFIN_EPS
                                and x[0] < xmax - DOLFIN_EPS
                                and x[1] > ymin + DOLFIN_EPS
                                and x[1] < ymax - DOLFIN_EPS
                                and (self.dim == 3
                                     and (x[2] > zmin + DOLFIN_EPS
                                          and x[2] < zmax - DOLFIN_EPS)))


# Outflow boundary
class OutflowBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], xmax)


class PsiMarker(Expression):

    def __init__(self, dim):
        self.dim = dim

    def eval(self, values, x):
        object = ObjectBoundary(self.dim)

        if(object.inside(x, True)):
            values[0] = 1.0
        else:
            values[0] = 0.0


# Problem definition
class Problem(ProblemBase):

    '''
        Provides the 2D/3D flow around a bluff body, where the bluff body is
        either a cylinder or a rectangular box. These problems were described in
        "Benchmark Computations of Laminar Flow Around a Cylinder" by M.
        Schaefer and S. Turek.
    '''

    def __init__(self, options, cube=False):

        global xmax, xcenter, ycenter

        ProblemBase.__init__(self, options)

        # Load mesh
        self.dim = int(options['dim'])
        self.Nx = options['Nx']
        options['Ny'] = None
        options['Nz'] = None
        self.cube = cube

        try:
            self.nu = options['nu']
        except:
            self.nu = 1E-3

        self.t0 = 0.
        try:
            self.T = options['T']
        except:
            self.T = 10

        # setup our domain and conditions
        if options['initial_mesh'] is not None:
            domain = options['initial_mesh']
            self.mesh = Mesh(domain)
        else:
            if self.dim == 2:  # 2D problem
                channel = Rectangle(Point(xmin, ymin), Point(xmax, ymax))
                if cube:
                    bluff = Rectangle(Point(xcenter - radius, ycenter - radius),
                                      Point(xcenter + radius, ycenter + radius))
                else:
                    bluff = Circle(Point(xcenter, ycenter), radius)
            else:  # 3D problem
                xmax = 2.5
                xcenter = 0.5
                channel = Box(Point(xmin, ymin, zmin), Point(xmax, ymax, zmax))
                if cube:
                    bluff = Box(Point(xcenter - radius, ycenter - radius, zmin),
                                Point(xcenter + radius, ycenter + radius, zmax))
                else:
                    bluff = Cylinder(Point(xcenter, ycenter, zmax),
                                     Point(xcenter, ycenter, zmin),
                                     radius, radius)

            domain = channel - bluff
            self.mesh = generate_mesh(domain, self.Nx)

        self.k = options['k']

        self.Re = self.Ubar*Diameter/self.nu

    def initial_conditions(self, W):
        if self.dim == 2:
            w0 = InitialConditions2D()
        else:
            w0 = InitialConditions3D()

        w0 = project(w0, W)

        return w0

    def boundary_conditions(self, W, t):
        self.U.t = t
        # Create inflow boundary condition
        subU = W.sub(0)
        bc0 = DirichletBC(subU, self.U, InflowBoundary())

        # Create no-slip boundary condition
        bc1 = DirichletBC(subU, self.noSlip,
                          NoSlipBoundary(self.dim))

        # Create outflow boundary condition for pressure
        # bc2 = DirichletBC(W.sub(1), Constant(0), OutflowBoundary())

        # Collect boundary conditions
        bcs = [bc0, bc1]

        return bcs

    def F(self, t):
        # forcing function for the momentum equation
        if self.dim == 2:
            f = Constant((0, 0))
        else:
            f = Constant((0, 0, 0))

        return f

    def update(self, W, t):
        # update the bc for each time step

        return self.boundary_conditions(W, t)

    def marker(self):
        # provide a marker to indicate if we are on the object
        return PsiMarker(self.dim)

    def functional(self, W, w):
        '''
            This is the functional used for adaptivity.
            We assume the problem is much like NSE.
        '''
        if W.mesh().topology().dim() == 2:
            (u, p, tau) = w.split()
        else:
            (u, p, tau) = w.split()

        # n = FacetNormal(mesh)
        # marker = problem.marker()

        M = u[0] * dx  # Mean of the x-velocity in the whole domain
        # M = marker*p*n[0]*ds  # Drag (only pressure)

        return M

    def __str__(self):
        return 'OldBCylinder'
