from pylab import *
import os

if len(sys.argv) < 3:
    print("Usage: python plot.py <adaptive data> <uniform data> title")
    print("Usage: python plot.py <adaptive data> title")
    exit(0)
elif len(sys.argv) < 4:
    titlestr = str(sys.argv[2])
    sys.argv[2] = None
else:
    titlestr = str(sys.argv[3])

adaptfile = sys.argv[1]
uniformfile = sys.argv[2]

rcParams['font.size'] = 16
rcParams['axes.titlesize'] = "large"
rcParams['axes.labelsize'] = "large"
rcParams['figure.figsize'] = 8, 5
rc('text', usetex=True)

adapt_data = loadtxt(adaptfile)
if uniformfile is not None:
    uniform_data = loadtxt(uniformfile)

adaptdofs = adapt_data[0:len(adapt_data) - 2, 0]
adapterr = abs(adapt_data[0:len(adapt_data) - 2, 1] - adapt_data[-1, 1])
adapterr_est = adapt_data[0:len(adapt_data) - 2, 2]

if uniformfile is not None:
    uniformdofs = uniform_data[0:len(uniform_data) - 2, 0]
    uniformerr = abs(uniform_data[0:len(uniform_data) - 2, 1]
                     - adapt_data[-1, 1])

hold(True)
title(titlestr, va='baseline', ha='center', loc='center')

plot(log10(adaptdofs), log10(adapterr), 'k')

if uniformfile is not None:
    plot(log10(uniformdofs), log10(uniformerr), 'r')

plot(log10(adaptdofs), log10(adapterr_est), 'k--')

if uniformfile is not None:
    legend(["Goal Oriented Error", "Uniform Error", "``True'' Error"],
           loc='best')
else:
    legend(["Goal Oriented Error", "``True'' Error"], loc='best')

ylabel('log10(error)')
xlabel('log10(dofs)')
grid(True)
savefig(os.path.splitext(os.path.basename(adaptfile))[0], bbox_inches='tight')
show()
