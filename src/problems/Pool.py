__author__ = "Robin Goix <robin.goix.rg@gmail.com>"
__date__ = "2014-04-22"
__license__  = "GNU GPL version 3 or any later version"
#
#   adapted from Square.py originally developed
#   by Erich L Foster <erichlf@gmail.com>
#
'''
    Provides the Pool Problem for the Peregrine System, where a moving object
    creates a disturbance producing a wave.
'''

from problembase import *
from numpy import array
from numpy import concatenate

x0 = -6.
x1 = 60.
y0 = -25.
y1 = 25.

#start location of our object
objectLeft = 0. #left bound
objectRight = 4. #right bound
objectBottom = -4. #lower bound
objectTop = -objectBottom #upper bound
channelBottom = objectBottom - 0.5
channelTop = objectTop + 0.5

#Physical Parameters
hd = 2. #Depth of the pool in the central lane [m]
hb = 0.3 #Depth at the boundaries [m]
ad = -0.3 #height of the moving object [m]
g = 9.8 #Gravity
lambda0 = 20. #typical wavelength
a0 = 0.8 #Typical wave height
h0 = 2. #Typical depth

#scaling our domain
x0 = x0/lambda0
x1 = x1/lambda0
y0 = y0/lambda0
y1 = y1/lambda0
objectLeft = objectLeft/lambda0
objectRight = objectRight/lambda0
objectBottom = objectBottom/lambda0
objectTop = objectTop/lambda0
channelBottom = channelBottom/lambda0
channelTop = channelTop/lambda0

#Scaled Parameters
hd = hd/h0 #depth
ad = ad/a0 #height of the moving object
hb = hb/h0 #depth at the boundary
c0 = (h0*g)**(0.5)

#maximum velocity of wave object
vmax = ((hb+a0)*g)**(0.5) #Max Speed of the moving object [m.s^(-1)]

class InitialConditions(Expression):
    def eval(self,values,x):
        values[0] = 0.
        values[1] = 0.
        values[2] = 0.

    def value_shape(self):
      return (3,)

# Slip boundary
class Y_SlipBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and \
               (x[1] < y0 + DOLFIN_EPS or x[1] > y1 - DOLFIN_EPS)

class X_SlipBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and \
               (x[0] < x0 + DOLFIN_EPS or x[0] > x1 - DOLFIN_EPS)

class Object(Expression):
    '''
        This expression takes the CST parameters and outputs the y cooresponding to
        x for the wave object. We use 3rd order Bernstein polynomials and thus we
        require 4 w-parameter, while the N parameters define the type of shape.
        The following code has been adapted from the code
        written by Pramudita Satria Palar for matlab which can be found at
        http://www.mathworks.com/matlabcentral/fileexchange/42239-airfoil-generation-using-cst-parameterization-method
    '''
    def __init__(self, vmax, H, N1, N2, W1, W2, W3, W4, dz, t):
        self.H = H #Object height
        self.N1 = N1 #shape parameter for 'back' of object
        self.N2 = N2 #shape parameter for 'front' of object
        self.w = (W1, W2, W3, W4) #weights for arc between 'front' and 'back'
        self.dz = dz #'leading' edge thickness
        self.t = t
        self.vmax = vmax

    def eval(self, value, x):
        H = self.H
        N1 = self.N1
        N2 = self.N2
        w = self.w
        dz = self.dz/H
        t = self.t
        vmax = self.vmax

        u = vmax*tanh(t)

        X = ((x[0] - objectLeft - u*t)/(objectRight - objectLeft), \
                (x[1] - objectBottom)/(objectTop - objectBottom))

        if(X[0]>=0 and X[0]<=1 and X[1]>=0 and X[1]<=1):
            # Class function; taking input of N1 and N2
            C = X[0]**N1*(1. - X[0])**N2

            # Shape function; using Bernstein Polynomials
            n = len(w) - 1 # Order of Bernstein polynomials

            S = 0;
            for i in range(0,n):
                K = float(factorial(n)/(factorial(i)*factorial(n-i)));
                S = S + w[i]*K*X[0]**i*(1. - X[0])**(n-i)

            value[0] = H*(C*S + X[0]*dz)*(1 - x[1]**2/objectTop**2)/0.25
        else:
            value[0] = 0.

    def value_shape(self):
        return ()

# Problem definition
class Problem(ProblemBase):
    '''
        Provides the Pool Problem for the Peregrine System, where a moving object
        creates a disturbance producing a wave.
    '''

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        #Scaling Parameters
        self.sigma = h0/lambda0
        self.epsilon = a0/h0

        # Create mesh
        self.Nx = options["Nx"]
        self.Ny = options["Ny"]
        self.mesh = RectangleMesh(x0, y0, x1, y1, self.Nx, self.Ny)

        self.mesh = self.Refine(self.mesh)

        #Scaled Parameters
        self.t0 = 0.
        self.T = options['T']*c0/lambda0 #Final time
        self.k = options['dt']*c0/lambda0 #time step

        #set up the CST shape parameterization
        ObjH = Constant(ad, name='ObjHeight')
        N1 = Constant(1., name='N1')
        N2 = Constant(1., name='N2')
        W1 = Constant(1., name='W1')
        W2 = Constant(1., name='W2')
        W3 = Constant(1., name='W3')
        W4 = Constant(1., name='W4')
        dz = Constant(0., name='dz')

        self.params = ['ObjHeight', 'N1', 'N2', 'W1', 'W2', 'W3', 'W4', 'dz']

        self.zeta0 = Object(vmax=vmax, H=ObjH, N1=N1, N2=N2, W1=W1, W2=W2, W3=W3, W4=W4, \
                dz=dz, t=self.t0)

        #Defintion of the shape of the seabed
        M1 = (hb - hd)/(y1 - objectTop)
        M2 = (hb - hd)/(y0 - objectBottom)
        M = 'x[1] > ' + str(objectTop) + ' ? ' + str(M1) + \
            ' : (x[1] < ' +  str(objectBottom) + ' ? ' + str(M2) + \
            ' : 0) '
        X = 'x[1] > ' + str(objectTop) + ' ? ' + str(objectTop) + \
            ' : (x[1] < ' +  str(objectBottom) + ' ? ' + str(objectBottom) + \
            ' : 0) '
        self.D = str(hd) + '+ (' + M + ')*(x[1] - (' + X + '))'

    def update_bathymetry(self, Q, t):
        '''
            Tells the solver what the bathymetry looks like at time t.
        '''
        D = Expression(self.D, element=Q.ufl_element(), annotate=False)

        self.zeta0.t = t
        zeta = project(self.zeta0, Q, annotate=False)
        self.zeta0.t = max(t - self.k, self.t0)
        zeta_ = project(self.zeta0, Q, annotate=False)
        self.zeta0.t = max(t - 2*self.k, self.t0)
        zeta__ = project(self.zeta0, Q, annotate=False)

        H = project(D + self.epsilon*zeta, Q, annotate=False)
        H_ = project(D + self.epsilon*zeta_, Q, annotate=False)

        return H, H_, zeta, zeta_, zeta__

    def Refine(self, mesh):
        '''
            Refine the mesh along the object's trajectory
        '''
        for i in range(0,4):
            cell_markers = CellFunction("bool", mesh)
            cell_markers.set_all(False)
            for cell in cells(mesh):
               p = cell.midpoint()
               if(p.y() > objectBottom*(1. - 1./(2.**(i+1)*lambda0)) and p.y() < objectTop*(1. + 1./(2.**(i+1)*lambda0))):
                    cell_markers[cell] = True
            mesh = refine(mesh, cell_markers)

        return mesh

    def initial_conditions(self, W):
        w0 = InitialConditions()
        w0 = project(w0,W,annotate=False)

        return w0

    def boundary_conditions(self, W, t):
        V = W.sub(0)

        # Create no-slip boundary condition for velocity
        bc_X = DirichletBC(V.sub(0), 0.0, X_SlipBoundary())
        bc_Y = DirichletBC(V.sub(1), 0.0, Y_SlipBoundary())

        bcs = [bc_X, bc_Y]

        return bcs

    def F1(self, t):
        #forcing function for the momentum equation
        return Expression(self.options['F1'],t=t)

    def F2(self, t):
        #mass source for the continuity equation
        return Expression(self.options['F2'],t=t)

    def __str__(self):
        return 'Pool'
