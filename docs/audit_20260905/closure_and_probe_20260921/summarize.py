#!/usr/bin/env python3
"""Reduce one closure_block.txt to the tables of the deliverable.
Reads only what the run printed; computes nothing the run did not."""
import re, sys, collections

path = sys.argv[1]
lines = [l.rstrip('\n') for l in open(path)]

def nums(s):
    return [float(x) for x in re.findall(r'[-+]?\d\.\d+E[-+]\d+', s)]

kk = collections.OrderedDict()   # k -> dict
act = {}                         # (k,h,dir) -> dict
order_h, order_dir, order_k = [], [], []

for l in lines:
    if ' K k=' in l:
        k = int(re.search(r'K k=(-?\d+)', l).group(1))
        if k not in order_k: order_k.append(k)
        d = kk.setdefault(k, {})
        d['completed'] = int(re.search(r'completed (\d+)', l).group(1))
        d['admissible'] = re.search(r'admissible (\w)', l).group(1)
        d['roots'] = int(re.search(r'roots_missing (\d+)', l).group(1))
        d['refusal'] = int(re.search(r'refusal (\d+)', l).group(1))
        m = re.search(r'increment\s+([-+]?\d\.\d+E[-+]\d+)', l)
        d['increment'] = m.group(1) if m else 'not taken'
    elif ' DIFF k=' in l:
        k = int(re.search(r'DIFF k=(-?\d+)', l).group(1))
        kk.setdefault(k, {})['diff_r'] = nums(l)
    elif ' DIFFJ k=' in l:
        k = int(re.search(r'DIFFJ k=(-?\d+)', l).group(1))
        kk.setdefault(k, {})['diff_j'] = nums(l)
    elif ' DIFFMAX k=' in l:
        k = int(re.search(r'DIFFMAX k=(-?\d+)', l).group(1))
        v = nums(l); cells = [int(x) for x in re.findall(r'E[-+]\d+\s+(\d+)', l)]
        kk.setdefault(k, {})['diffmax'] = (v, cells)
    elif ' CLOSURE k=' in l:
        k = int(re.search(r'CLOSURE k=(-?\d+)', l).group(1))
        kk.setdefault(k, {})['closure'] = nums(l)[0]
        kk[k]['channel'] = l.split('channel ')[-1]
    elif ' CHANNELS k=' in l:
        k = int(re.search(r'CHANNELS k=(-?\d+)', l).group(1))
        kk.setdefault(k, {})['channels'] = nums(l)
    else:
        m = re.match(r' \(closure_probe\) (\w+) k=(-?\d+) h=\s*([-\d.E+]+) dir=(\S+)(.*)', l)
        if not m: continue
        tag, k, h, dr, rest = m.group(1), int(m.group(2)), float(m.group(3)), m.group(4), m.group(5)
        if h not in order_h: order_h.append(h)
        if dr not in order_dir: order_dir.append(dr)
        e = act.setdefault((k, h, dr), {})
        if tag == 'ACT':
            e['ok'] = re.search(r'ok=(\w\w\w)', rest).group(1)
            v = nums(rest)
            e['eps_nom'], e['eps'], e['dispY'], e['dispD'] = v[0], v[1], v[2], v[3]
            sw = re.search(r'swept (\d+) (\d+)', rest)
            e['swept'] = (int(sw.group(1)), int(sw.group(2)))
            ro = re.search(r'roots (\d+) (\d+)', rest)
            e['roots'] = (int(ro.group(1)), int(ro.group(2)))
            bl = re.search(r'blocked (\d+) (\w)', rest)
            e['blocked'], e['backward'] = int(bl.group(1)), bl.group(2)
        elif tag in ('JPROD', 'JFWD', 'JCTR', 'JFC'):
            e[tag] = nums(rest)
        elif tag in ('DH', 'DK'):
            v = nums(rest)
            e[tag] = v[:5]; e[tag + 'rel'] = v[5]
        elif tag == 'JHASH':
            hs = re.findall(r'fwd ([0-9A-F]{16}) ctr ([0-9A-F]{16})', rest)
            if hs: e['hash'] = hs[0]

hs = sorted(order_h)
print('## The closure count, at the state\n')
print('| k | passes completed | increment the sweep exited on | admissible | roots missing | refusal | closure increment against the converged closure | channel | ||Dr^-1 (F-F_jac)||_2 mass / momentum / energy / all | ||Dj^-1 (F-F_jac)||_2 all |')
print('|---|---|---|---|---|---|---|---|---|---|')
for k in order_k:
    d = kk[k]
    kn = 'converged' if k == 0 else str(k)
    dr = d.get('diff_r', [0]*5); dj = d.get('diff_j', [0]*5)
    print('| %s | %d | %s | %s | %d | %d | %.3E | %s | %.4E / %.4E / %.4E / %.4E | %.4E |' %
          (kn, d['completed'], d['increment'], d['admissible'], d['roots'], d['refusal'],
           d.get('closure', 0.0), d.get('channel', ''), dr[0], dr[1], dr[2], dr[4], dj[4]))

for dr in order_dir:
    print('\n## The action along `%s`: the two dimensions\n' % dr)
    print('| k \\ h | ' + ' | '.join('%g' % h for h in hs) + ' |')
    print('|---|' + '---|'*len(hs))
    for k in order_k:
        kn = 'converged' if k == 0 else str(k)
        row = []
        for h in hs:
            e = act.get((k, h, dr))
            row.append('%.3E' % e['JCTR'][4] if e and 'JCTR' in e else '-')
        print('| %s | %s |' % (kn, ' | '.join(row)))
    print('\nThe distance from the action of the same k at the production arc h = 1, relative:\n')
    print('| k \\ h | ' + ' | '.join('%g' % h for h in hs) + ' |')
    print('|---|' + '---|'*len(hs))
    for k in order_k:
        kn = 'converged' if k == 0 else str(k)
        row = []
        for h in hs:
            e = act.get((k, h, dr))
            row.append('%.2E' % e['DHrel'] if e and 'DHrel' in e else '-')
        print('| %s | %s |' % (kn, ' | '.join(row)))
    print('\nThe distance from the action at the converged closure, same arc, relative:\n')
    print('| k \\ h | ' + ' | '.join('%g' % h for h in hs) + ' |')
    print('|---|' + '---|'*len(hs))
    for k in order_k:
        kn = 'converged' if k == 0 else str(k)
        row = []
        for h in hs:
            e = act.get((k, h, dr))
            row.append('%.2E' % e['DKrel'] if e and 'DKrel' in e else '-')
        print('| %s | %s |' % (kn, ' | '.join(row)))

print('\n## Endpoint record: passes completed at Y+eps v and Y-eps v, roots missing, blocked components\n')
print('| k | h | direction | ok (probe, +, -) | eps requested | eps used | displacement ||D^-1 eps v|| | passes (+, -) | roots (+, -) | blocked | backward |')
print('|---|---|---|---|---|---|---|---|---|---|---|')
for k in order_k:
    for h in hs:
        for dr in order_dir:
            e = act.get((k, h, dr))
            if not e or 'ok' not in e: continue
            kn = 'converged' if k == 0 else str(k)
            print('| %s | %g | %s | %s | %.3E | %.3E | %.3E | %d, %d | %d, %d | %d | %s |' %
                  (kn, h, dr, e['ok'], e['eps_nom'], e['eps'], e['dispD'],
                   e['swept'][0], e['swept'][1], e['roots'][0], e['roots'][1],
                   e['blocked'], e['backward']))
