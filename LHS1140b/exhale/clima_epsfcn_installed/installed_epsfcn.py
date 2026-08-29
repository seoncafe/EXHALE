"""Read the finite-difference step out of the installed Clima library.

`AdiabatClimate_simple_solver` passes `epsfcn` to MINPACK by reference, so the
value the installed build actually uses is a constant in `.rodata`.  This
disassembles the routine, takes the last three stack arguments pushed before
the `hybrd` call -- `ml`, `mu`, `epsfcn` in reverse push order -- and prints
the double the third of them points at, together with `factor`.

usage: python installed_epsfcn.py [path to _clima*.so]
"""
import os
import re
import struct
import subprocess
import sys

EX = '/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00'
SO = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    EX, 'env/photochem/lib/python3.11/site-packages/photochem',
    '_clima.cpython-311-x86_64-linux-gnu.so')

data = open(SO, 'rb').read()

sections = []
for line in subprocess.run(['readelf', '-S', '-W', SO],
                           capture_output=True, text=True).stdout.splitlines():
    if line.strip().startswith('['):
        q = line[line.index(']') + 1:].split()
        if len(q) >= 5:
            try:
                sections.append((q[0], int(q[2], 16), int(q[3], 16), int(q[4], 16)))
            except ValueError:
                pass


def double_at(va):
    for _, addr, off, size in sections:
        if addr and addr <= va < addr + size:
            fo = off + (va - addr)
            return struct.unpack('<d', data[fo:fo + 8])[0]
    return None


sym = subprocess.run(['nm', SO], capture_output=True, text=True).stdout
m = re.search(r'^([0-9a-f]+) t (\S*adiabatclimate_simple_solver\S*)$',
              sym, re.M)
if not m:
    sys.exit('AdiabatClimate_simple_solver not found in %s' % SO)
start = int(m.group(1), 16)
print('symbol   %s at 0x%x' % (m.group(2), start))

asm = subprocess.run(['objdump', '-d', '--start-address=0x%x' % start,
                      '--stop-address=0x%x' % (start + 0x1000), SO],
                     capture_output=True, text=True).stdout.splitlines()

call = None
pushes = []          # rip-relative target of each push, in code order
rax = None           # target of the most recent `lea X(%rip),%rax`, unconsumed
for i, line in enumerate(asm):
    if re.search(r'callq.*minpack_module_MOD_hybrd1>', line):
        print('WARNING: this build still calls hybrd1')
    if re.search(r'callq.*minpack_module_MOD_hybrd>', line):
        call = i
        break
    g = re.search(r'lea\s+-?0x[0-9a-f]+\(%rip\),%rax.*#\s([0-9a-f]+)', line)
    if g:
        rax = int(g.group(1), 16)
        continue
    if re.search(r'\bpush\s+%rax\b', line):
        pushes.append(rax)
        rax = None
    elif re.search(r'\bpush\b', line):
        pushes.append(None)

if call is None:
    sys.exit('no call to hybrd found: this build does not carry the repair')

# stack arguments 7, 8, 9 ... are pushed last first: the final push is arg 7.
# hybrd(fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode, factor, ...)
args = list(reversed(pushes))          # args[0] = arg 7 = ml
names = ['ml', 'mu', 'epsfcn', 'diag', 'mode', 'factor', 'nprint']
for k, name in enumerate(names):
    if k >= len(args):
        break
    va = args[k]
    if va is None:
        continue
    print('arg %-2d %-7s -> 0x%x  = %r' % (k + 7, name, va, double_at(va)))
