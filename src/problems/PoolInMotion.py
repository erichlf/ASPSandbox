__author__ = "Robin Goix <robin.goix.rg@gmail.com>"
__date__ = "2014-04-22"
__license__  = "GNU GPL version 3 or any later version"
#
#   adapted from Square.py originally developed 
#   by Erich L Foster <erichlf@gmail.com>
#

'''
This is a simple problem with no forcing and on a rectangle.
This is basically a blank problem that we can adapt with optional inputs.
'''

from problembase import *
from numpy import array

x0 = -30.
x1 = 10.
y0 = -15.
y1 = 15.
# No-slip boundary

class Y_SlipBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and \
               (x[1] < y0 + DOLFIN_EPS or x[1] > y1 - DOLFIN_EPS)
           
class VelocityStream_Boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and \
               (x[0] < x0 + DOLFIN_EPS or x[0] > x1 - DOLFIN_EPS)

class Entry_Boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (x[0] > x1 - DOLFIN_EPS)
           
# Problem definition
class Problem(ProblemBase):
#   2D channel flow.

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        global x0, x1, y0, y1

        # Create mesh
        Nx = options["Nx"]
        Ny = options["Ny"]
        lambda0 = float(options["lambda0"])
        x0 = x0/lambda0
        x1 = x1/lambda0
        y0 = y0/lambda0
        y1 = y1/lambda0
        mesh = RectangleMesh(x0, y0, x1, y1, Nx, Ny)
        
        #Refine the mesh along the object's trajectory
        cell_markers = CellFunction("bool", mesh)
        cell_markers.set_all(False)

        for cell in cells(mesh):
            p = cell.midpoint()
            if p.y() > -3.5/lambda0 and p.y() < 3.5/lambda0: 
                cell_markers[cell] = True
            
        self.mesh = refine(mesh, cell_markers)
        
        cell_markers2 = CellFunction("bool", self.mesh)
        cell_markers2.set_all(False)

        for cell in cells(self.mesh):
            p = cell.midpoint()
            if p.y() > -3./lambda0 and p.y() < 3./lambda0 and p.x() > -12./lambda0:
                cell_markers2[cell] = True
            
        self.mesh = refine(self.mesh, cell_markers2)
        
        cell_markers3 = CellFunction("bool", self.mesh)
        cell_markers3.set_all(False)
        
        for cell in cells(self.mesh):
            p = cell.midpoint()
            if p.y() > -3./lambda0 and p.y() < 3./lambda0 and p.x() > -10./lambda0 and p.x() < 6./lambda0:
                cell_markers3[cell] = True
            
        self.mesh = refine(self.mesh, cell_markers3)
        
        cell_markers4 = CellFunction("bool", self.mesh)
        cell_markers4.set_all(False)
        
        for cell in cells(self.mesh):
            p = cell.midpoint()
            if p.y() > -2.5/lambda0 and p.y() < 2.5/lambda0 and p.x() > -8./lambda0 and p.x() < 4./lambda0:
                cell_markers4[cell] = True
            
        self.mesh = refine(self.mesh, cell_markers4)
           
        g = 9.8 #Gravity
        lambda0 = self.options['lambda0'] #typical wavelength
        a0 = self.options['a0'] #Typical wave height
        h0 = self.options['h0'] #Typical depth
        c0 = (h0*g)**(0.5)
        U = (g*(h0+0.4))**0.5
        self.U = h0/a0*U/c0
        
        
    def initial_conditions(self, V, Q):
        u0 = Expression(("-U", "0.0"),U=self.U)
        eta0 = Expression("0.0")

        return u0, eta0

    def boundary_conditions(self, V, Q, t):
        # Create no-slip boundary condition for velocity
        bc_X_u = DirichletBC(V, Expression(("-U","0.0"),U=self.U), VelocityStream_Boundary())
        bc_X_eta = DirichletBC(Q, 0.0, Entry_Boundary())
        bc_Y_u = DirichletBC(V.sub(1), 0.0, Y_SlipBoundary())
        
        bcs = [bc_X_u, bc_Y_u, bc_X_eta]
        
        return bcs

    def F1(self, t):
        #forcing function for the momentum equation
        return Expression(self.options['F1'],t=t)

    def F2(self, t):
        #mass source for the continuity equation
        return Expression(self.options['F2'],t=t)

    def __str__(self):
        return 'PoolInMotion'
