__author__ = "Kristian Valen-Sendstad <kvs@simula.no>"
__date__ = "2009-10-01"
__copyright__ = "Copyright (C) 2009-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg, 2010.
# Modified by Erich L Foster, 2013

from problembase import *
from numpy import array

# Constants related to the geometry
bmarg = 1.e-3 + DOLFIN_EPS
xmin = 0.0; xmax = 2.2
ymin = 0.0; ymax = 0.41
zmin = 0.0; zmax = 0.41
xcenter = 0.2; ycenter = 0.2
radius = 0.05; Diameter = 2.*radius
Um = 1.5 #max velocity

class InitialConditions2D(Expression):
    def eval(self, values, x):
        values[0] = 0.
        values[1] = 0.
        values[2] = 0.

    def value_shape(self):
        return (3,)

class InitialConditions3D(Expression):
    def eval(self,values,x):
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
    def __init__(self, dim, cube):
        SubDomain.__init__(self)
        self.dim = dim
        self.cube = cube

    def inside(self, x, on_boundary):
        dx = x[0] - xcenter
        dy = x[1] - ycenter
        r = sqrt(dx*dx + dy*dy)

        Cube = near(x[0], xcenter - radius) or near(x[0], xcenter + radius) \
                or near(x[1], ycenter - radius) or near(x[1], ycenter + radius)

        return on_boundary \
                and (near(x[1], ymin) or near(x[1], ymax) \
                or (self.dim == 3 and (near(x[2], zmin) or near(x[2], zmax))) \
                or (not self.cube and r < radius + bmarg) \
                or (self.cube and Cube))

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
        ProblemBase.__init__(self, options)

        global xmax, xcenter, ycenter

        # Load mesh
        self.dim = int(options['dim'])
        self.Nx = options['Nx']
        options['Ny'] = None
        options['Nz'] = None
        self.cube = cube

        #setup our domain and conditions
        if self.dim == 2: #2D problem
            channel = Rectangle(xmin, ymin, xmax, ymax)
            if cube:
                bluff = Rectangle(xcenter - radius, ycenter - radius, \
                        xcenter + radius, ycenter + radius)
            else:
                bluff = Circle(xcenter, ycenter, radius)
            self.noSlip = Constant((0,0))
            self.U = Expression(('4*Um*x[1]*(H - x[1])/(H*H)', '0.0'), Um=Um, H=ymax)
            Ubar = 2.*self.U((0,ymax/2.))[0]/3.
        else: #3D problem
            xmax = 2.5
            xcenter = 0.5;
            channel = Box(xmin, ymin, zmin, xmax, ymax, zmax)
            if cube:
                bluff = Box(xcenter - radius, ycenter - radius, zmin, \
                        xcenter + radius, ycenter + radius, zmax)
            else:
                bluff = Cylinder(Point(xcenter, ycenter, zmax), \
                        Point(xcenter, ycenter, zmin), radius)
            self.noSlip = Constant((0,0,0))
            self.U = Expression(('16*Um*x[1]*x[2]*(H - x[1])*(H - x[2])/pow(H,4)', \
                    '0.0', '0.0'), Um=Um**2, H=ymax)
            Ubar = 4.*self.U((0,ymax/2.,zmax/2.))[0]/9.

        domain = channel - bluff
        self.mesh = Mesh(domain, self.Nx)

        #rescale Reynolds number to the problem
        options['Re'] = Ubar*Diameter*options['Re']

        self.t0 = 0.
        self.T = options['T']
        self.k = options['dt']

        #since Cube relies on this code we need a lot of selfs
        self.channel = channel

    def initial_conditions(self, W):
        if self.dim == 2:
            w0 = InitialConditions2D()
        else:
            w0 = InitialConditions3D()

        w0 = project(w0,W)

        return w0

    def boundary_conditions(self, W, t):
        # Create inflow boundary condition
        bc0 = DirichletBC(W.sub(0), self.U, InflowBoundary())

        # Create no-slip boundary condition
        bc1 = DirichletBC(W.sub(0), self.noSlip, NoSlipBoundary(self.dim, self.cube))

        # Create outflow boundary condition for pressure
        bc2 = DirichletBC(W.sub(1), Constant(0), OutflowBoundary())

        # Collect boundary conditions
        bcs = [bc0, bc1, bc2]

        return bcs

    def F1(self, t):
        #forcing function for the momentum equation
        if self.dim == 2:
            f = Constant((0,0))
        else:
            f = Constant((0,0,0))

        return f

    def F2(self, t):
        #mass source for the continuity equation
        return Constant(0)

    def __str__(self):
        return 'Cylinder'
