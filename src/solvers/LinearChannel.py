#!/usr/bin/python

from dolfin import *
import sys

set_log_active(False)

#define parameters
def params():
    h = 1.63 #depth
    g = 9.81 #gravity
    f0 = 6.16E-5 #reference Coriolis paramter
    beta = 2.07E-11 #beta plane approximation
    return h, g, f0, beta


def f(f0,beta): #Coriolis parameter
    return Expression('f0 + beta*x[1]', f0=f0, beta=beta)

#define domain 
def domain():
    x0 = -1
    y0 = -1
    x1 = 1
    y1 = 1
    return x0, x1, y0, y1

#define mesh size
Nx = 64 #number of elements in x-direction 
Ny = 64 #number of elements in y-direction

#define order of mesh
P_u = 1 #order of the velocity elements
P_eta = 1 #order of the pressure elements

#define time interval, time step, and solver 
t = 0 #start time 
T = 1 #final time
dt = 0.01 #time step
theta = 0.5 #time stepping method

#define the pressure boundary conditions
p_in = Expression('sin(pi*t)',t=t)
p_out = 0.0

#where to save solution
directory = 'LinearChannel'
hBaseName = 'Height'
UBaseName = 'Velocity'
ext = 'pvd'
hFile = File(directory + '/' + hBaseName + '.' + ext)
UFile = File(directory + '/' + UBaseName + '.' + ext)

#create mesh and define function space
x0, x1, y0, y1 = domain()
mesh = RectangleMesh(x0,y0,x1,y1,Nx,Ny)
#mesh = Mesh('NorthAtlanticV2.xml')

# Class representing the initial conditions
class InitialConditions(Expression):
    def eval(self, values, x):
        h, g, f0, beta = params()
        x0, x1, y0, y1 = domain()
        a = 0.0
        mx = (x1 + x0)/2.
        my = (y1 + y0)/2.
        sx = 5.92E-2
        sy = 5.92E-2
        values[0] = 0.0#2*g*a*b*x[1]*exp(-b*(x[0]**2 + x[1]**2))
        values[1] = 0.0#-2*g*a*b*x[0]*exp(-b*(x[0]**2 + x[1]**2))
        values[2] = a*exp(-((x[0]-mx)**2/(2*sx**2) + (x[1]-my)**2/(2*sy**2)))
    def value_shape(self):
        return (3,)

def noslip(x, on_boundary):
    x0, x1, y0, y1 = domain()
    return on_boundary

def inflow(x, on_boundary):
    x0, x1, y0, y1 = domain()
    return on_boundary and near(x[0], x0)

def outflow(x, on_boundary):
    x0, x1, y0, y1 = domain()
    return on_boundary and near(x[0], x1)

# Class for interfacing with the Newton solver
class LinearSWE(NonlinearProblem):
    def __init__(self, w, w_k, v, chi, t, bcs):
        NonlinearProblem.__init__(self)
        h, g, f0, beta = params()
        U, eta = split(w)
        U_k, eta_k = split(w_k)

        # eta_(k+theta)
        eta_mid = (1.0-theta)*eta_k + theta*eta

        # U_(k+theta)
        U_mid = (1.0-theta)*U_k + theta*U

        #weak form of the equations
        L0 = (eta - eta_k)*chi*dx \
            - dt*h*inner(U_mid,grad(chi))*dx 
        L1 = inner(U - U_k,v)*dx \
            + f(f0,beta)*dt*(U_mid[0]*v[1] - U_mid[1]*v[0])*dx \
            + dt*g*inner(grad(eta_mid),v)*dx 
        L = L0 + L1

        # Compute directional derivative about w in the direction of dw (Jacobian)
        a = derivative(L, w, dw)

        self.L = L
        self.a = a
        self.bcs = bcs
        self.reset_sparsity = True
    def F(self, b, x):
        assemble(self.L, tensor=b, bcs=self.bcs)
    def J(self, A, x):
        assemble(self.a, tensor=A, reset_sparsity=self.reset_sparsity,
            bcs=self.bcs)
        self.reset_sparsity = False

#compiler options
parameters["form_compiler"]["optimize"]     = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "quadrature"

#define our functional spaces
V = VectorFunctionSpace(mesh, 'CG', P_u)
Q = FunctionSpace(mesh, 'CG', P_eta)
W = MixedFunctionSpace([V, Q])

#define trial and test function
dw = TrialFunction(W) #direction of the Gateaux derivative
(v, chi) = TestFunctions(W)

#boundary conditions
no_slip  = DirichletBC(V, (0, 0), noslip)
in_flow = DirichletBC(Q, p_in, inflow)
out_flow = DirichletBC(Q, p_out, outflow)

bcs = [no_slip, in_flow, out_flow]

#define functions 
w = Function(W) #solution from current step
w_k = Function(W) #solution from previous step

U = Function(V)
U_k = Function(V)
eta = Function(Q)
eta_k = Function(Q)

#split mixed functions
(du, deta) = split(dw) #direction for the Gateaux derivative

#setup initial conditions
w_init = InitialConditions()
w = project(w_init,W)
w_k = project(w_init,W)

#save the initial condition
hFile << w.split()[1] #height
UFile << w.split()[0] #height

# Create nonlinear problem and Newton solver
solver = NewtonSolver("lu")
prm = solver.parameters
prm['convergence_criterion'] = 'incremental' 
prm['absolute_tolerance'] = 1E-6
prm['relative_tolerance'] = 1E-5
prm['maximum_iterations'] = 10
prm['relaxation_parameter'] = 1.0
prm['report'] = False

#initialize plot
#viz = plot(w_k.split()[0], mesh=mesh, 
#    title='Height')

# Compute solution
while t < T:
    t += dt

    w_k.vector()[:] = w.vector() #set the solution for the previous time step

    if(t <= 1.0): #update boundary pressure
        p_in.t = t

    problem = LinearSWE(w, w_k, v, chi, t, bcs) #build the Shallow Water Equations FE

    sys.stdout.flush()
    sys.stdout.write("t=%G\r" % t)
    solver.solve(problem, w.vector()) #solve our problem

    hFile << w.split()[1] #height
    UFile << w.split()[0] #height

    # Plot solution and mesh
#    viz.plot(w.split()[0]) 
