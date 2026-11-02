"""Second pass on the (Delta,w) grid: refit each point of grid_<part>.json from an extra start
(the w=20 keV solution, which lies in a different basin with chi_N+2<0) and keep the lower chi2."""
import sys, json
sys.path.insert(0, '.')
import cmfit as C
part = int(sys.argv[1])
pts = json.load(open(f'grid_{part}.json'))
alt = json.load(open('fit_w20.json'))['P']
FREE = ['eN', 'kN', 'cN', 'kI', 'cI', 'qN', 'qI']
out = []
for p in pts:
    if p['w'] < 12:
        out.append(dict(p, polished=False)); continue
    st = {k: alt[k] for k in FREE}
    P, c2, err = C.fit(FREE, fixed=dict(D=p['D'], w0=p['w'], w2=p['w']), starts=[st], max_nfev=60, tag=f"polish D={p['D']} w={p['w']}")
    if c2 < p['chi2']: out.append(dict(D=p['D'], w=p['w'], chi2=c2, P=P, polished=True))
    else: out.append(dict(p, polished=False))
    json.dump(out + pts[len(out):], open(f'grid_{part}p.json', 'w'))
json.dump(out, open(f'grid_{part}p.json', 'w'))
print('DONE polish', part)
