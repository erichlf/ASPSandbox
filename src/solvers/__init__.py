__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2018-12-26"
__license__ = "GNU GPL version 3 or any later version"

# List of solvers
solvers = [
    'NSE',
    'SWE',
    'Peregrine',
    'QGE',
    'QGECGDG',
    'DensityNSE',
    'Heat',
    'ADR',
    'MixedWaveALE',
    'OldroydB',
    'ThreeField_Stokes',
    'Poisson',
    'SS_NSE',
    'SS_3F_Stokes',
]


def Solver(name, options):  # Wrapper for solver classes
    # Return solver instance for given solver name"
    exec('from .{} import Solver as NamedSolver'.format(name), globals())
    return NamedSolver(options)
