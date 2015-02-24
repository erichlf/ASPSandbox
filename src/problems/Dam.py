__author__ = "Johan Jansson <jjan@csc.kth.se>"
__date__ = "2013-08-27"
__license__ = "GNU GPL version 3 or any later version"

from AFES import *
from AFES import Problem as ProblemBase
from mshr import *

geps = 1.0e-8

d = 1.
eta = 0.1
g = 9.8

xMin, xMax = -2 * d, 2 * d
yMin, yMax = -d, d

xMinDam, xMaxDam = -0.2 * d, 0.2 * d
yMinDam, yMaxDam = -d, -0.6 * d

rhoMin, rhoMax = 1., 1000.
At = (rhoMax - rhoMin) / (rhoMin + rhoMax)  # Atwood number


class InitialConditions(Expression):

    def eval(self, values, x):
        values[0] = 0.
        values[1] = 0.
        if x[1] <= -0.7 * d or (x[0] >= -1.8 * d and x[0] <= -d
                                and x[1] <= -0.2 * d):
            values[2] = rhoMax
        else:
            values[2] = rhoMin

        values[3] = 0.

    def value_shape(self):
        return (4,)


class NoSlipBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary


class DamBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary and (x[0] <= 2 * d - geps
                                and x[0] >= -2 * d + geps
                                and x[1] <= d - geps and x[1] >= - d + geps)


class Problem(ProblemBase):

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        self.mesh = self.Mesh(options)

        self.t0 = 0.
        try:
            self.T = options['T']
        except:
            self.T = 10

        try:
            self.nu = option['nu']
        except:
            self.nu = 1E-3

        self.At = At
        self.Ubar = 1.
        try:
            self.CFL = options['CFL']
        except:
            self.CFL = 5.
        self.k = self.time_step(self.Ubar, self.mesh)  # mesh size

        # what functional should we use?
        fDic = {'pressure_drag': self.pDrag,
                'density': self.density,
                }
        try:
            self.functional = fDic[options['functional']]
        except:
            self.functional = fDic['pressure_drag']

    def Mesh(self, options):  # setup our domain
        self.Nx = options['Nx']
        if options['initial_mesh'] is not None:
            mesh = Mesh(initial_mesh)
        else:
            channel = Rectangle(Point(xMin, yMin), Point(xMax, yMax))
            dam = Rectangle(Point(xMinDam, yMinDam), Point(xMaxDam, yMaxDam))

            domain = channel - dam
            mesh = generate_mesh(domain, self.Nx)

        return mesh

    def time_step(self, u, mesh):
        return self.CFL * mesh.hmin()/u

    def initial_conditions(self, W):

        self.artificial_viscosity(W)

        # wf = File("wb.pvd")
        # wf << self.wb
        w0 = InitialConditions()
        # w0 = project(w0,W)
        ww0 = Function(W)
        wt = TestFunction(W)
        F = (inner(ww0, wt)
             + 1E-3 * inner(grad(ww0), grad(wt))
             - inner(w0, wt))*dx
        solve(F == 0, ww0)

        return ww0

    def artificial_viscosity(self, W):
        V = FunctionSpace(W.mesh(), 'CG', 1)

        bc0 = DirichletBC(V, Constant(1.0), NoSlipBoundary())

        h = CellSize(V.mesh())

        wb = Function(V)
        wt = TestFunction(V)

        F = wb*wt*dx + h*inner(grad(wb), grad(wt))*dx

        solve(F == 0, wb, bc0)
        self.wb = wb

    def boundary_conditions(self, W, t):
        # Create no-slip boundary condition for velocity
        bc0 = DirichletBC(W.sub(0), Constant((0.0, 0.0)), NoSlipBoundary())
        # Create zero boundary condition for pressure
        # bc1 = DirichletBC(W.sub(2), Constant(0.0), NoSlipBoundary())
        bcs = [bc0]  # , bc1]

        return bcs

    def F(self, t):
        # forcing function for the momentum equation

        return Expression(('0.0', '-g'), g=g, t=t)

    def measure(self, W):  # create a subdomain for integrating over the circle
        dam = DamBoundary()
        boundaries = FacetFunction("size_t", W.mesh())
        boundaries.set_all(0)
        dam.mark(boundaries, 1)

        return Measure("ds")[boundaries]

    def pDrag(self, W, w):  # drag on dam object

        p = w[3]
        n = FacetNormal(W.mesh())

        # obtain a measure for the various domains ds(1) is the dam
        ds = self.measure(W)

        return p * n[0] * ds(1)  # Drag (only pressure)

    def density(self, W, w):
        rho = w[2]

        beta = Expression(('50.0*exp(-50.0*(pow(x[0] - 1.0, 2) '
                           + ' + pow(x[1] - -0.5, 2)))'))

        return rho * beta * dx

    def __str__(self):
        return 'Dam'
