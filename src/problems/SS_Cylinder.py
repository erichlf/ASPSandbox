from Cylinder import *
from Cylinder import Problem as Cylinder


class Problem(Cylinder):

    '''
        Purely a container to call Cylinder. This
        convert Cylinder.py to a steady state problem.
    '''

    def __init__(self, options, cube=False):
        Cylinder.__init__(self, options, cube=cube)

        H = ymax
        if self.mesh.topology().dim() == 2:
            self.noSlip = Constant((0, 0))
            self.U = Expression(('4*Um*x[1]*(H - x[1])/(H*H)', '0.0'),
                                Um=Um, H=ymax)
            self.Ubar = 4. / 3. * Um * ymax * (H - ymax / 2.) / (H * H)
        else:
            self.noSlip = Constant((0, 0, 0))
            self.U = Expression(('16*Um*x[1]*x[2]*(H - x[1])*(H - x[2])' +
                                 '/pow(H,4)', '0.0', '0.0'),
                                Um=Um ** 2, H=ymax)
            self.Ubar = 16. / 9. * Um * ymax * zmax * (H - ymax / 2.) \
                * (H - zmax / 2.) / pow(H, 4)

        self.Re = self.Ubar*Diameter/self.nu

    def boundary_conditions(self, W):
        # Create inflow boundary condition
        bc0 = DirichletBC(W.sub(0), self.U, InflowBoundary())

        # Create no-slip boundary condition
        bc1 = DirichletBC(W.sub(0), self.noSlip, NoSlipBoundary(self.dim))

        # Create outflow boundary condition for pressure
        bc2 = DirichletBC(W.sub(1), Constant(0), OutflowBoundary())

        # Collect boundary conditions
        bcs = [bc0, bc1, bc2]

        return bcs

    def F(self):
        return Constant((0.0, 0.0))

    def __str__(self):
        return 'SS_Cylinder'
