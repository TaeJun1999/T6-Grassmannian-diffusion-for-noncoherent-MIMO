#!/usr/bin/env python3
"""Q-12 (b) part 3: a TRUE NLOS point (verified with a LOS-only solve) on Munich, clutter-free and with car clutter."""
import json, numpy as np, time
import t6_q12_rt_cars as C
from t6_q12_metrics import metrics
st=(56.0,16.0)
# probe direction (specular only, no cars)
best=None
for dd in [(1,0),(-1,0),(0,1),(0,-1)]:
    pws=[np.linalg.norm(C.H_at((st[0]+dd[0]*s, st[1]+dd[1]*s), False)[0])**2 for s in np.linspace(0,3,4)]
    sp=10*np.log10(max(pws)/max(min(pws),1e-30)) if min(pws)>0 else 99
    if best is None or sp<best[1]: best=(dd,sp)
d=best[0]; print('NLOS_true direction',d,round(best[1],1),flush=True)
out={'scene':'munich, true NLOS (56,16)','BS':C.BS,'runs':[]}
def run_nocars(name,Delta,B):
    t0=time.time(); H=[]; los=[]
    for b in range(B):
        Hb,n,l=C.H_at((st[0]+d[0]*C.V*b*Delta, st[1]+d[1]*C.V*b*Delta)); H.append(Hb); los.append(l)
    H=np.stack(H); m=metrics(H,W=16,lags=(1,2,4,8,16,32,47)); Hs,_,_=C.H_at(st,False)
    m.update({'name':name,'start':list(st),'direction':list(d),'Delta_s':Delta,'B':B,'clutter':'none, S=0.4','dist_per_block_lambda':C.V*Delta/C.LAM,
              'los_fraction_blocks':float(np.mean(los)),'specular_power_fraction':float(np.linalg.norm(Hs)**2/max(np.linalg.norm(H[0])**2,1e-30)),
              'mean_power':float(np.mean(np.linalg.norm(H,axis=(1,2))**2)),'time_s':time.time()-t0,'H_re':H.real.tolist(),'H_im':H.imag.tolist()})
    print(f"{name} no-cars D={Delta*1e3:g}ms: los={m['los_fraction_blocks']:.2f} spec_frac={m['specular_power_fraction']:.2f} rhoG(1)={m['rho_G']['1']:.3f} rhoG(4)={m['rho_G']['4']:.3f} rhoS(1)={m['rho_S']['1']:.3f} r95={m['r_eff95']} PEF2(1)={m['PEF_r2'].get('1',float('nan')):.3f} PEF2(32)={m['PEF_r2'].get('32',float('nan')):.3f} PEFe(32)={m['PEF_reff95'].get('32',float('nan')):.3f} LT(32)={m['LT_overlap_reff95'].get('32',float('nan')):.3f} power={m['mean_power']:.2e} ({m['time_s']:.0f}s)",flush=True)
    return m
for Delta,B in [(10e-3,80),(2e-3,80)]:
    out['runs'].append(run_nocars('NLOS_true',Delta,B)); json.dump(out,open('t6_q12_rt_nlos.json','w'))
C.place_cars(st,d)
for Delta,B in [(10e-3,80),(2e-3,80)]:
    m=C.run('NLOS_true',st,d,Delta,B); m['mean_power']=float(np.mean(np.linalg.norm(np.array(m['H_re'])+1j*np.array(m['H_im']),axis=(1,2))**2)); out['runs'].append(m); json.dump(out,open('t6_q12_rt_nlos.json','w'))
print('done')
