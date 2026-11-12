"""Recompute every number and regenerate every generated figure of the two
papers.  Usage:  python run_all.py            (checks + figures, ~2 minutes)
                 python run_all.py --checks   (numbers only)
                 python run_all.py --figures  (figures only)
Exit code 0 iff every check passes."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, 'checks'))


def main(argv):
    do_checks = '--figures' not in argv
    do_figs = '--checks' not in argv
    ok = True
    if do_checks:
        import numcheck
        ok = numcheck.run(verbose=True)
        print('NUMCHECK ALL PASS' if ok else 'NUMCHECK: SOME FAIL')
    if do_figs:
        sys.path.insert(0, os.path.join(HERE, 'figures'))
        import shadow_and_chord, prefix_and_switching, windows_and_column, dens_smoothed, dual_step4
        shadow_and_chord.main()
        prefix_and_switching.main()
        windows_and_column.main()
        dens_smoothed.main()
        dual_step4.main()
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
