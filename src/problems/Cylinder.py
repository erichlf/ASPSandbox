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
import sys

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

    def value_shape(self):
        return (3,)


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

        b = on_boundary and (x[0] > xmin + DOLFIN_EPS
                             and x[0] < xmax - DOLFIN_EPS
                             and x[1] > ymin + DOLFIN_EPS
                             and x[1] < ymax - DOLFIN_EPS)
        if self.dim == 3:
            b = b and x[2] > zmin + DOLFIN_EPS and x[2] < zmax - DOLFIN_EPS

        return b


# Outflow boundary
class OutflowBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], xmax)


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

        H = ymax
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

        try:
            varying = options['varying']
        except:
            varying = False

        if self.mesh.topology().dim() == 2:
            self.noSlip = Constant((0, 0))
            if varying:
                self.U = Expression(('4*Um*x[1]*(H - x[1])*sin(pi*t/8)/(H*H)',
                                     '0.0'), Um=Um, H=ymax, t=self.t0)
                self.T = 8
            else:
                self.U = Expression(('4*Um*x[1]*(H - x[1])/(H*H)', '0.0'),
                                    Um=Um, H=ymax, t=self.t0)
            self.Ubar = 4. / 3. * Um * ymax * (H - ymax / 2.) / (H * H)
        else:
            self.noSlip = Constant((0, 0, 0))
            self.U = Expression(('16*Um*x[1]*x[2]*(H - x[1])*(H - x[2])' +
                                 '/pow(H,4)', '0.0', '0.0'),
                                Um=Um ** 2, H=ymax, t=self.t0)
            self.Ubar = 16. / 9. * Um * ymax * zmax * (H - ymax / 2.) \
                * (H - zmax / 2.) / pow(H, 4)

        try:
            self.CFL = options['CFL']
        except:
            self.CFL = 100.
        self.k = self.time_step(self.Ubar, self.mesh)  # mesh size

        self.Re = self.Ubar*Diameter/self.nu

        # what functional should we use?
        fDic = {'velocity': self.velocity,
                'drag': self.drag,
                'pressure_drag': self.pDrag,
                'lift': self.lift,
                'pressure_lift': self.pLift,
                }
        self.CdMax = 0
        try:
            self.functional = fDic[options['functional']]
        except:
            self.functional = fDic['velocity']

    def initial_conditions(self, W):
        if self.dim == 2:
            w0 = InitialConditions2D()
        else:
            w0 = InitialConditions3D()

        w0 = project(w0, W)

        return w0

    def time_step(self, u, mesh):
        return self.CFL * mesh.hmin()/u

    def boundary_conditions(self, W, t):

        self.U.t = t
        # Create inflow boundary condition
        bc0 = DirichletBC(W.sub(0), self.U, InflowBoundary())

        # Create no-slip boundary condition
        bc1 = DirichletBC(W.sub(0), self.noSlip, NoSlipBoundary(self.dim))

        # Create outflow boundary condition for pressure
        bc2 = DirichletBC(W.sub(1), Constant(0), OutflowBoundary())

        # Collect boundary conditions
        bcs = [bc0, bc1, bc2]

        return bcs

    def F1(self, t):  # forcing function for the momentum equation
        if self.dim == 2:
            f = Constant((0, 0))
        else:
            f = Constant((0, 0, 0))

        return f

    def F2(self, t):  # mass source for the continuity equation
        return Constant(0)

    def update(self, W, t):  # update the bc for each time step

        return self.boundary_conditions(W, t)

    def velocity(self, W, w):  # Mean of the x-velocity in the whole domain
        '''
            This is the functional used for adaptivity.
            We assume the problem is much like NSE.
        '''
        if W.mesh().topology().dim() == 2:
            u = as_vector((w[0], w[1]))
        else:
            u = as_vector((w[0], w[1], w[2]))

        return u[0] * dx

    def measure(self, W):  # create a subdomain for integrating over the circle
        bluff = ObjectBoundary(self.dim)
        boundaries = FacetFunction("size_t", W.mesh())
        boundaries.set_all(0)
        bluff.mark(boundaries, 1)

        return Measure("ds")[boundaries]

    def bodyForce(self, W, w, lift=False):  # Calculate the force on the bluff body
        if W.mesh().topology().dim() == 2:
            (u, p) = (as_vector((w[0], w[1])), w[2])
        else:
            (u, p) = (as_vector((w[0], w[1], w[2])), w[3])

        n = FacetNormal(W.mesh())
        I = Identity(2)
        sigma = p*I - self.nu*self.epsilon(u)
        if lift:
            theta = Constant((0.0, 1.0))
        else:
            theta = Constant((1.0, 0.0))

        ds = self.measure(W)

        return dot(dot(sigma, n), theta) * ds(1)

    def drag(self, W, w):  # drag functional
        '''
            This is the functional used for adaptivity.
            We assume the problem is much like NSE.
        '''
        if W.mesh().topology().dim() == 2:
            return 2. / (self.Ubar**2. * Diameter) * self.bodyForce(W, w)
        else:
            return 2. / (self.Ubar**2. * Diameter * H) * self.bodyForce(W, w)

    def pDrag(self, W, w):  # functional for Drag (only pressure)
        '''
            This is the functional used for adaptivity.
            We assume the problem is much like NSE.
        '''
        if W.mesh().topology().dim() == 2:
            p = w[2]
        else:
            p = w[3]

        n = FacetNormal(W.mesh())

        ds = self.measure(W)

        return p*n[0]*ds(1)

    def lift(self, W, w):  # drag functional
        '''
            This is the functional used for adaptivity.
            We assume the problem is much like NSE.
        '''
        if W.mesh().topology().dim() == 2:
            return 2. / (self.Ubar**2. * Diameter) \
                * self.bodyForce(W, w, lift=True)
        else:
            return 2. / (self.Ubar**2. * Diameter * H) \
                * self.bodyForce(W, w, lift=True)

    def pLift(self, W, w):  # functional for Lift (only pressure)
        '''
            This is the functional used for adaptivity.
            We assume the problem is much like NSE.
        '''
        if W.mesh().topology().dim() == 2:
            p = w[2]
        else:
            p = w[3]

        n = FacetNormal(W.mesh())

        ds = self.measure(W)

        return p * n[1] * ds(1)  # Lift (only pressure)

    def epsilon(self, z):
        return 0.5*(grad(z) + grad(z).T)

    def __str__(self):
        return 'Cylinder'
