__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__ = "GNU GPL version 3 or any later version"

# List of problems
problems = [
    'Bubble',
    'Cavity',
    'Channel',
    'Cube',
    'Cylinder',
    'Dam',
    'DamBreak',
    'DoubleGyre',
    'Drop',
    'Duct',
    'Ground',
    'Hole',
    'Point',
    'Pool',
    'SmallHole',
    'Square',
    'TwoFluid',
    'SS_Cylinder',
    'OldBCylinder',
    'SS_Hole',
    'SS_3F_Cylinder',
]


def Problem(name, options):  # Wrapper for problem classes
    # Return problem instance for given problem name"
    exec('from %s import Problem as NamedProblem' % name)
    return NamedProblem(options)
