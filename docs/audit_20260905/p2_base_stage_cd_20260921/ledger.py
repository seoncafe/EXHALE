import re, sys
fn=sys.argv[1]
solve=0; rec={}; out=[]
pat_it=re.compile(r'\(JFNK\) it\s+(\d+)\s+\|\|R\|\|=\s*(\S+)\s+\|\|Fs\|\|2=\s*(\S+)\s+lam=\s*(\S+)\s+dtau=\s*(\S+)\s+gm=\s*(\d+)\s+worst r=\s*(\S+)\s+worst row: (.*?)\s+solve=(\d+)')
for line in open(fn):
    m=re.search(r'forcing term\s*(\S+): the cycle reached\s*(\S+) in (\d+) product\(s\) of (\d+); (.*)', line)
    if m: rec['eta'],rec['gmres'],rec['np'],rec['m'],rec['why']=m.groups(); continue
    m=re.search(r'rows entering the step \(mass, momentum, energy\):\s+(\S+)\s+(\S+)\s+(\S+)', line)
    if m: rec['in']=m.groups(); continue
    m=re.search(r'linear residual left, over the row of F \(mass, momentum, energy\):\s+(\S+)\s+(\S+)\s+(\S+)', line)
    if m: rec['rlin']=m.groups(); continue
    m=re.search(r'rows now \(mass, momentum, energy\):\s+(\S+)\s+(\S+)\s+(\S+)', line)
    if m: rec['now']=m.groups(); continue
    m=re.search(r'judged slots .*? of the state:\s+(\S+)\s+(\S+)\s+(\S+)\s+of the trial:\s+(\S+)\s+(\S+)\s+(\S+)', line)
    if m: rec['js']=m.groups(); continue
    m=re.search(r'cells outside the tolerance of their row, of 500: mass (\d+) \(worst cell (\d+)\), momentum (\d+) \((\d+)\), energy (\d+) \((\d+)\)', line)
    if m: rec['out']=m.groups(); continue
    m=pat_it.search(line)
    if m:
        it,rn,f2,lam,dtau,gm,wr,wrow,sv=m.groups()
        out.append(dict(solve=int(sv),it=int(it),rnorm=rn,f2=f2,lam=lam,dtau=dtau,gm=gm,worst=wrow.strip(),**rec))
        rec={}
print("solve it | eta  gmres_rel np | rows_in m,p,e | r_lin/F m,p,e | rows_now m,p,e | judged state->trial | lam dtau ||R|| ||Fs||2 | worst | cells_out m,p,e")
for d in out:
    print("%d %3d | %s %s %s | %s | %s | %s | %s -> %s | %s %s %s %s | %s | %s"%(
        d['solve'],d['it'],d.get('eta','-'),d.get('gmres','-'),d.get('np','-'),
        ' '.join(d.get('in',('-','-','-'))),' '.join(d.get('rlin',('-','-','-'))),
        ' '.join(d.get('now',('-','-','-'))),
        d.get('js',('-',)*6)[0],d.get('js',('-',)*6)[3],
        d['lam'],d['dtau'],d['rnorm'],d['f2'],d['worst'],
        ','.join(d.get('out',('-',)*6)[i] for i in (0,2,4))))
