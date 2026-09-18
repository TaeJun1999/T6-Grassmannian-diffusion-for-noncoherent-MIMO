#!/usr/bin/env python3
"""Table of all Sionna-RT trajectory runs (t6_q12_rt.json, t6_q12_rt_cars.json, t6_q12_rt_nlos.json) with corrected point labels."""
import json, os, numpy as np
LAB = {'LOS_A': 'P1 (40,48) LOS 49m, LOS 88%', 'NLOS_A': 'P2 (-8,48) LOS 41m, LOS 54%+refl 27/14%', 'NLOS_B': 'P3 (8,40) LOS 32m, LOS 84%',
       'LOS_far': 'P4 (72,72) LOS 85m, LOS 83%', 'NLOS_true': 'P5 (56,16) NLOS 48m (verified), 4 spec paths'}
rows = []
for f in ['t6_q12_rt_metrics.json', 't6_q12_rt.json', 't6_q12_rt_cars.json', 't6_q12_rt_nlos.json']:
    if not os.path.exists(f): continue
    d = json.load(open(f))
    if f != 't6_q12_rt_metrics.json' and os.path.exists('t6_q12_rt_metrics.json'): continue  # metrics file supersedes
    for r in d['runs']:
        rows.append((LAB[r['name']], r.get('clutter', 'none'), r.get('S', 0.4), r['Delta_s'] * 1e3, r['dist_per_block_lambda'], r['specular_power_fraction'],
                     r['rho_G']['1'], r['rho_G']['4'], r['rho_G'].get('16', float('nan')), r['rho_S']['1'], r['r_eff95'],
                     r['PEF_r2'].get('1', float('nan')), r['PEF_r2'].get('32', float('nan')), r['PEF_reff95'].get('32', float('nan')), r['LT_overlap_reff95'].get('32', float('nan')),
                     10 * np.log10(r.get('mean_power', np.nan) / 1e-9) if r.get('mean_power') else float('nan')))
L = []
L.append("Sionna RT 2.1.0, scene 'munich', BS (8.5,21,27) ULA N=32, UE ULA M=2 @1.5 m, 28 GHz, v=3 m/s, max_depth 3, los+specular+diffuse, seed 1")
L.append("NOTE: the earlier 'NLOS_A/NLOS_B' labels were wrong (LOS test used normalized delays); P1-P4 are LOS points, P5 is verified NLOS.")
L.append("power = mean ||H||_F^2 in dB relative to 1e-9 (array sum, unit-gain iso elements)")
L.append(f"{'point':<44} {'clutter':<26} {'S':>3} {'D[ms]':>5} {'d/lam':>5} {'spec%':>5} {'rG(1)':>6} {'rG(4)':>6} {'rG(16)':>6} {'rS(1)':>6} {'r95':>3} {'PEF2(1)':>7} {'PEF2(32)':>8} {'PEFe(32)':>8} {'LT(32)':>6} {'pow':>6}")
for r in rows:
    L.append(f"{r[0]:<44} {r[1]:<26} {r[2]:>3} {r[3]:>5g} {r[4]:>5.2f} {100*min(r[5],1):>5.0f} {r[6]:>6.3f} {r[7]:>6.3f} {r[8]:>6.3f} {r[9]:>6.3f} {r[10]:>3d} {r[11]:>7.3f} {r[12]:>8.3f} {r[13]:>8.3f} {r[14]:>6.3f} {r[15]:>6.1f}")
L.append("")
L.append("Decision rule (open_questions Q-12, fixed 2026-09-17 before these runs): gain decorrelated <=> rho_G(1) <= 0.5 for v*Delta/lambda in [0.3, 3];")
L.append("LT persistence <=> LT_overlap(32) >= 0.8 and PEF_reff(32) >= 0.85.  Verdict per run:")
for r in rows:
    gd = (0.3 <= r[4] <= 3) and r[6] <= 0.5; lt = (r[14] >= 0.8) and (r[13] >= 0.85)
    L.append(f"  {r[0][:28]:<28} {r[1][:18]:<18} D={r[3]:g}ms: gain-decorrelated={'YES' if gd else 'no '} ({r[6]:.2f})  LT-persistent={'YES' if lt else 'no '} (LT {r[14]:.2f}, PEFe {r[13]:.2f})")
open('t6_q12_rt_summary.txt', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
