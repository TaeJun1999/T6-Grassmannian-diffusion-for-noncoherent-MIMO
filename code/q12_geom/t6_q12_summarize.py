#!/usr/bin/env python3
"""Summaries for Q-12 (a): t6_q12_geom.json, t6_q12_longterm.json, t6_ceiling_LT.json (+ project t6_ceiling_results.json)."""
import json, sys, numpy as np
L=[]
def P(s=''): L.append(s); print(s)
g=json.load(open('t6_q12_geom.json')); lt=json.load(open('t6_q12_longterm.json'))
P("=== Q-12(a) part 1: inter-block gain correlation vs instantaneous receive-subspace overlap (t6_q12_geom.json, seed %d, B=%d, runs=%d) ==="%(g['seed'],g['B'],g['n_runs']))
P("motion label: phi_v=90 -> radial (UE at broadside moving away, AoA fixed); phi_v=0 -> tangential (AoA drifts at v/R)")
for name,c in g['cases'].items():
    kw=c['kwargs']; P(f"-- {name}  {kw}")
    P(f"{'Delta[ms]':>9} {'d/lam':>6} {'rhoG(1)':>8} {'|J0|':>6} {'rhoS(1)':>8} {'rhoS(8)':>8} {'rhoS(16)':>9} {'rhoS(32)':>9} {'top2frac':>9}")
    for D,r in c['results'].items():
        P(f"{float(D)*1e3:9.2f} {r['dist_per_block_lambda']:6.2f} {r['rho_G']['1']:8.3f} {r['jakes_ref']['1']:6.3f} {r['rho_S']['1']:8.3f} {r['rho_S']['8']:8.3f} {r['rho_S']['16']:9.3f} {r['rho_S']['32']:9.3f} {r['top2_energy_fraction_window']:9.3f}")
P()
P("=== Q-12(a) part 2: long-term subspace (t6_q12_longterm.json; W=16-block windows; PEF = predictive energy fraction from the PAST window) ===")
P(f"{'case':<42} {'r90':>4} {'r95':>4} {'PR':>5} {'rhoS1':>6} {'PEF2(1)':>8} {'PEF2(16)':>9} {'PEF2(48)':>9} {'PEFe(1)':>8} {'PEFe(48)':>9} {'LT(16)':>7} {'LT(48)':>7}")
for name,c in lt['cases'].items():
    P(f"{name:<42} {c['r_eff90']:4.1f} {c['r_eff95']:4.1f} {c['participation_ratio']:5.2f} {c['rho_S_inst_lag1']:6.3f} {c['PEF_r2']['1']:8.3f} {c['PEF_r2']['16']:9.3f} {c['PEF_r2']['48']:9.3f} {c['PEF_reff95']['1']:8.3f} {c['PEF_reff95']['48']:9.3f} {c['LT_overlap_reff95']['16']:7.3f} {c['LT_overlap_reff95']['48']:7.3f}")
P()
d=json.load(open('t6_ceiling_LT.json')); base=json.load(open('/mnt/project/t6_ceiling_results.json'))
pf={c['N']:c for c in base['configs']}
P("=== Revised gate-(0) ceiling with a LONG-TERM-subspace genie (t6_ceiling_LT.json, seed 20260917; reduced problem (T,2,r_eff,2) at SNR + 10log10(N/r_eff)) ===")
P("SER 1e-2 operating points [dB]; 'LT-genie' = struct-ML of the reduced problem shifted by -10log10(N/r_eff) (isotropic 2-plane inside the r_eff-dim LT subspace = pessimistic genie)")
P(f"{'regime':<10} {'r_eff':>5} {'K':>5} {'prior-free struct':>18} {'LT-genie':>9} {'exact-U genie':>14} {'LT gap':>7} {'exact gap':>10}")
for c in d['configs']:
    T,N=c['T'],c['N']; Nf=32 if T==16 else 16; sh=10*np.log10(Nf/N)
    for K in ('64','1024'):
        r=c['finite'][K+'_snr_at_ser1e-2']; f=pf[Nf]['finite'][K+'_snr_at_ser1e-2']
        P(f"{'(%d,2,%d,2)'%(T,Nf):<10} {N:5d} {K:>5} {f['struct']:18.2f} {r['struct']-sh:9.2f} {f['genie']:14.2f} {f['struct']-(r['struct']-sh):7.2f} {f['struct']-f['genie']:10.2f}")
P()
P("Rate domain: not evaluated for the LT genie. The finite-K struct-ML MI saturates (K cap) at the main-regime SNRs, and the reduced problem's")
P("continuous-USTM quantities (R_genie' = exact-U genie, GMI_svd' = lower bound) do not bracket the LT-genie rate tightly; a K=8192 low-SNR run")
P("of the reduced problem would be needed (Q-08/Q-14). Judgment unit for (0) stays dB @ SER 1e-2 (D-07).")
open('t6_q12_summary.txt','w').write('\n'.join(L)+'\n')
