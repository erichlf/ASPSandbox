__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2015-08-10"
__license__ = "GNU GPL version 3 or any later version"

from ASP import *
from ASP import Problem as ProblemBase
from mshr import *

d = 0.38  # height of water column below dam in cm
d0 = 15.  # height of the water column behind dam in cm
g = 9.8  # gravity

# channel definition
xMin, xMax = 0., 220.
yMin, yMax = 0., 20.
xLock = 38.  # location of the lock gate in cm

rhoMin, rhoMax = 1., 1000.  # approximation of air and water
At = (rhoMax - rhoMin) / (rhoMin + rhoMax)  # Atwood number


class InitialConditions(Expression):

    def eval(self, values, x):
        values[0], values[1], values[3] = 0., 0., 0.
        if (x[0] <= 38 and x[1] <= 15) or (x[0] > 38 and x[1] <= 0.38):
            values[2] = rhoMax
        else:
            values[2] = rhoMin

    def value_shape(self):
        return (4,)


class Problem(ProblemBase):

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        self.mesh = self.Mesh(options)

        self.t0 = 0.
        try:
            self.T = options['T']
        except:
            self.T = 2

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
        self.functional = self.density

    def Mesh(self, options):  # setup our domain
        self.Nx = options['Nx']
        self.Ny = options['Ny']
        if options['initial_mesh'] is not None:
            mesh = Mesh(initial_mesh)
        else:
            mesh = RectangleMesh(xMin, yMin, xMax, yMax, self.Nx, self.Ny)

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

        bc0 = DirichletBC(V, Constant(1.0), 'on_boundary')

        h = CellSize(V.mesh())

        wb = Function(V)
        wt = TestFunction(V)

        F = (wb * wt + h * inner(grad(wb), grad(wt))) * dx

        solve(F == 0, wb, bc0)
        self.wb = wb

    def boundary_conditions(self, W, t):
        # Create no-slip boundary condition for velocity
        bc0 = DirichletBC(W.sub(0), Constant((0.0, 0.0)), 'on_boundary')
        bcs = [bc0]

        return bcs

    def F(self, t):
        # forcing function for the momentum equation

        return Expression(('0.0', '-g'), g=g, t=t)

    def density(self, W, w):  # density functional
        rho = w[2]

        beta = Expression(('50.0*exp(-50.0*(pow(x[0] - xLock, 2) '
                           + ' + pow(x[1] - d0, 2)))'), d0=d0, xLock=xLock)

        return rho * beta * dx

    def __str__(self):
        return 'DamBreak'
