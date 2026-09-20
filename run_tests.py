#!/usr/bin/env python
"""
Run the complete KEYLINK test suite:   python run_tests.py [-v] [-k pattern]

Equivalent to  python -m unittest discover -s tests -t .  (or  pytest tests  if pytest is installed).
The tests run on a two-year climate fixture (tests/fixtures/) and take about one minute.
Exit status is 0 when everything passes, 1 otherwise.
"""
import argparse
import os
import sys
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-v', '--verbose', action='store_true', help='list every test')
    ap.add_argument('-k', metavar='PATTERN', help='only run test files/classes/methods whose name contains PATTERN')
    args = ap.parse_args()

    os.chdir(HERE)
    sys.path.insert(0, HERE)
    os.environ.setdefault('MPLBACKEND', 'Agg')

    loader = unittest.TestLoader()
    if args.k:
        loader.testNamePatterns = ['*%s*' % args.k]
    suite = loader.discover(start_dir='tests', pattern='test_*.py', top_level_dir=HERE)
    n = suite.countTestCases()
    print('KEYLINK test suite: %d tests' % n)
    t0 = time.time()
    result = unittest.TextTestRunner(verbosity=2 if args.verbose else 1).run(suite)
    print('%s in %.1f s' % ('ALL TESTS PASSED' if result.wasSuccessful() else 'FAILURES', time.time()-t0))
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    sys.exit(main())
