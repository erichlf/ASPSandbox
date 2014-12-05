import sys
from time import time

from AFES import *
import Heat as Solver
import Square as Problem


def main(args):
    # Get options
    options = OPTIONS.copy()

    # Set debug level
    set_log_active(options["debug"])

    print 'Solving %s for the %s problem.' % (solver_name, problem_name)
    # Create problem and solver
    problem = Problem(options)
    solver = Solver(options)

    # Solve problem with solver
    wct = time()
    w = solver.solve(problem)

    # Compute elapsed time
    wct = time() - wct

    print 'Solved %s for the %s problem in %g seconds.' % (solver.__str__,
                                                           problem.__str__, wct)

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
