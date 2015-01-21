from AFES import *
from AFES import Problem as ProblemBase

ymax = 0.41
xmax = 2.2
xmin = 0.0
xmax = 2.2
ymin = 0.0
ymax = 0.41
xcenter = 0.2
ycenter = 0.2
radius = 0.05
Diameter = 2. * radius

'''
    SubDomains for defining boundary conditions
'''


class WallBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (near(x[1], ymax) or near(x[1], ymin))


class InnerBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and not (x[0] > xmin + DOLFIN_EPS
                                    or x[1] > ymin + DOLFIN_EPS
                                    or x[0] < xmax - DOLFIN_EPS
                                    or x[1] < ymax - DOLFIN_EPS)


class InflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], xmin)


class OutflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], xmax)


class Inflow(Expression):  # Coefficients for defining boundary conditions
    def eval(self, values, x):
        values[0] = 4*1.0*(x[1]*(ymax-x[1]))/(ymax*ymax)
        values[1] = 0

    def value_shape(self):
        return (2,)


class PsiMarker(Expression):
    def eval(self, values, x):
        ib = InnerBoundary()

        if(ib.inside(x, True)):
            values[0] = 1.0
        else:
            values[0] = 0.0


class Problem(ProblemBase):
    def __init__(self, options):
        ProblemBase.__init__(self, options)

        self.mesh = Mesh("data/turek_cylinder_01.xml")

        try:
            self.nu = options['nu']
        except:
            self.nu = 1E-3

    def boundary_conditions(self, W):
        u0_0 = Constant((0.0, 0.0))
        u0_0p = Constant(0.0)

        # primal solver boundary conditions
        bc_0_0 = DirichletBC(W.sub(0), u0_0, WallBoundary())
        bc_0_1 = DirichletBC(W.sub(0), u0_0, InnerBoundary())
        bc_0_2 = DirichletBC(W.sub(0), Inflow(), InflowBoundary())
        bc_0_3p = DirichletBC(W.sub(1), u0_0p, OutflowBoundary())

        return [bc_0_0, bc_0_1, bc_0_2, bc_0_3p]

    def F(self):
        return Constant((0.0, 0.0))

    def functional(self, W, w):

        nu = self.nu
        mesh = W.mesh()
        psimarker = PsiMarker()

        (u, p) = (as_vector((w[0], w[1])), w[2])

        n = FacetNormal(mesh)

        I = Identity(2)
        sigma = p*I - nu*self.epsilon(u)
        theta = Constant((1.0, 0.0))

        g = Expression(
            ("200.0*exp(-200.0*(pow(x[0] - 0.5, 2) + pow(x[1] - 0.3, 2)))",
             "0.0"))

        M1 = psimarker*p*n[0]*ds  # Drag (only pressure)
        M2 = psimarker*p*n[1]*ds  # Lift (only pressure)
        M3 = inner(g, u)*dx  # Mean of the velocity in a region
        M4 = psimarker*dot(dot(sigma, n), theta)*ds  # Drag (full stress)
        M5 = u[0]*dx  # Mean of the x-velocity in the whole domain

        return M5

    def epsilon(self, z):
        return 0.5*(grad(z) + grad(z).T)
