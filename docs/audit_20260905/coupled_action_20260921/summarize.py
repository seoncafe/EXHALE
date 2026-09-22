#!/usr/bin/env python3
"""Reduce one closure_probe block to the tables of the coupled Stage B record."""
import re, sys, collections

def load(path):
    rows = collections.defaultdict(dict)
    meta = []
    kt = {}
    for ln in open(path):
        if '(closure_probe)' not in ln: continue
        t = ln.split('(closure_probe)',1)[1].rstrip('\n')
        f = t.split()
        tag = f[0]
        if tag in ('K','CLOSURE','DIFF','DIFFMAX','DIFFJ','CHANNELS'):
            k = int(f[1].split('=')[1]); kt.setdefault(k,{})[tag] = t
            continue
        if tag in ('ACT','JPROD','JFWD','JCTR','JFC','DH','DK','JHASH','MOVED'):
            m = re.match(r'\s*\S+\s+k=(-?\d+)\s+h=\s*(\S+)\s+dir=(\S+)', t)
            k = int(m.group(1)); h = float(m.group(2)); d = m.group(3)
            rows[(k,h,d)][tag] = f
            continue
        meta.append(t)
    return meta, kt, rows

def num(f, i): return float(f[i])

meta, kt, rows = load(sys.argv[1])
print('---- identity ----')
for m in meta: print(m)
print()
print('---- the closure count at the state ----')
for k in sorted(kt):
    for tag in ('K','CLOSURE','DIFF','DIFFJ'):
        if tag in kt[k]: print(kt[k][tag])
print()
ks = sorted(set(k for (k,h,d) in rows))
hs = sorted(set(h for (k,h,d) in rows))
ds = []
for (k,h,d) in rows:
    if d not in ds: ds.append(d)
def cls5(f, tag):
    # ... ' ||...||_2 by class ' then 5 numbers
    return [float(x) for x in f[-5:]]
print('---- action (central) ||Dr^-1 J v||_2, whole vector ----')
for d in ds:
    for k in ks:
        line = ' '.join('%11.4e'%cls5(rows[(k,h,d)]['JCTR'],'')[4] if (k,h,d) in rows and 'JCTR' in rows[(k,h,d)] else '   -' for h in hs)
        print('%-12s k=%-2d %s'%(d,k,line))
print('h list: '+' '.join('%.0e'%h for h in hs))
print()
print('---- plateau in h: ||J(h)-J(1)|| / ||J(1)|| ----')
for d in ds:
    for k in ks:
        line = ' '.join('%11.3e'%float(rows[(k,h,d)]['DH'][-1]) if (k,h,d) in rows and 'DH' in rows[(k,h,d)] else '   -' for h in hs)
        print('%-12s k=%-2d %s'%(d,k,line))
print()
print('---- movement with k: ||J(k)-J(conv)|| / ||J(conv)|| ----')
for d in ds:
    for k in ks:
        line = ' '.join('%11.3e'%float(rows[(k,h,d)]['DK'][-1]) if (k,h,d) in rows and 'DK' in rows[(k,h,d)] else '   -' for h in hs)
        print('%-12s k=%-2d %s'%(d,k,line))
print()
print('---- forward against central over the rows the direction MOVES (mean / worst / n) ----')
for d in ds:
    for k in ks:
        cells=[]
        for h in hs:
            r = rows.get((k,h,d),{}).get('MOVED')
            if r is None: cells.append('       -'); continue
            # ... 'moved rows/mean/worst' n mean worst 'worst moved row' ...
            i = r.index('rows/mean/worst', r.index('moved'))
            cells.append('%9.3e/%-4s'%(float(r[i+2]), r[i+1]))
        print('%-12s k=%-2d %s'%(d,k,' '.join(cells)))
print()
print('---- admissibility and the endpoint record (ACT) ----')
bad=0
for key in sorted(rows):
    r = rows[key].get('ACT')
    if r is None: continue
    ok = r[r.index('ok=TTT') if 'ok=TTT' in r else 4]
    txt=' '.join(r)
    if 'ok=TTT' not in txt or 'roots 0 0' not in txt or 'blocked 0' not in txt:
        print('NOT CLEAN: '+txt); bad+=1
print('points with a non-clean endpoint record: %d of %d'%(bad,len([1 for k in rows if 'ACT' in rows[k]])))
