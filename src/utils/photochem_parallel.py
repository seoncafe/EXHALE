#!/usr/bin/env python3
"""Run the lower-atmosphere adapter with the existing OpenMP thread budget.

Photochem and clima limit OpenMP to one thread on import. Load those modules
before applying OMP_NUM_THREADS so their import defaults do not override the
requested budget. This changes execution of their existing parallel regions,
not the chemical reactions, radiative transfer, or convergence criteria.
"""
import os
from pathlib import Path
import runpy
import sys


def requested_threads():
    value = os.environ.get('OMP_NUM_THREADS', '1')
    try:
        count = int(value)
    except ValueError as exc:
        raise ValueError('OMP_NUM_THREADS must be one positive integer') from exc
    if count < 1 or str(count) != value.strip():
        raise ValueError('OMP_NUM_THREADS must be one positive integer')
    return count


def main():
    count = requested_threads()
    import photochem  # noqa: F401
    if '--climate' in sys.argv:
        import photochem.clima  # noqa: F401
    from threadpoolctl import threadpool_info, threadpool_limits
    with threadpool_limits(limits=count, user_api='openmp'):
        libraries = [row for row in threadpool_info() if row['user_api'] == 'openmp']
        if not libraries or any(row['num_threads'] != count for row in libraries):
            raise RuntimeError('The installed OpenMP libraries did not accept the Photochem thread budget')
        print('Photochem OpenMP threads: requested %d; active %s'
              % (count, ','.join(str(row['num_threads']) for row in libraries)), flush=True)
        adapter = Path(__file__).resolve().with_name('photochem_to_lower_profile.py')
        runpy.run_path(str(adapter), run_name='__main__')


if __name__ == '__main__':
    main()
