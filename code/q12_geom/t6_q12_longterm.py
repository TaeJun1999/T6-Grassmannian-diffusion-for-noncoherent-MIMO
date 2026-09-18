#!/usr/bin/env python3
"""
t6_q12_longterm.py -- Q-12 (a), part 2: instantaneous vs long-term receive subspace.
For each case and block spacing Delta: blocks b = 0..B-1, window W.
  r_eff(W): number of eigenvalues of C_W = sum_{b in window} H_b^H H_b needed for 90% / 95% energy; participation ratio.
  PEF_r(k): predictive energy fraction = E ||H_b P||_F^2 / ||H_b||_F^2, P = projector on the top-r eigenspace of the
            covariance of the PAST window [b-k-W, b-k)  (what a subspace tracker / long-term prior could supply).
            r = 2 (instantaneous-rank object of D-01) and r = r_eff95.
  LT overlap(k): overlap of top-r_eff subspaces of windows [0,W) and [k,k+W) (drift of the long-term subspace).
Seed 20260917. Imports Geometry from t6_q12_geom.py.
"""
import json, sys, time
import numpy as np
from t6_q12_geom import Geometry, SEED

def topk_proj(C, r):
    ev, V = np.linalg.eigh(C)
    return V[:, -r:]

def run(name, kw, Delta, B=80, W=16, n_runs=30, lags=(1, 4, 16, 32, 48)):
    rng = np.random.default_rng(SEED + hash(name) % 1000)
    out = {'kwargs': kw, 'Delta_s': Delta, 'W': W}
    acc = {'r90': 0, 'r95': 0, 'pr': 0.0, 'rhoS1': 0.0, 'top2_inst_in_LT': 0.0}
    pef2 = {k: [0.0, 0.0] for k in lags}; pefe = {k: [0.0, 0.0] for k in lags}; ltov = {k: [0.0, 0] for k in lags}
    reff_list = []
    for run_i in range(n_runs):
        geo = Geometry(rng, **kw)
        H = np.stack([geo.channel(b * Delta) for b in range(B)])
        H /= np.linalg.norm(H, axis=(1, 2), keepdims=True)
        C_all = np.einsum('bmn,bmk->nk', H.conj(), H)
        ev = np.sort(np.linalg.eigvalsh(C_all))[::-1]; cs = np.cumsum(ev) / ev.sum()
        r90 = int(np.searchsorted(cs, 0.90) + 1); r95 = int(np.searchsorted(cs, 0.95) + 1)
        acc['r90'] += r90; acc['r95'] += r95; acc['pr'] += ev.sum() ** 2 / (ev ** 2).sum()
        reff_list.append(r95)
        # instantaneous lag-1 overlap
        U = []
        for b in range(B):
            _, _, Vh = np.linalg.svd(H[b], full_matrices=False); U.append(Vh[:2].conj().T)
        acc['rhoS1'] += np.mean([np.linalg.norm(U[b].conj().T @ U[b + 1]) ** 2 / 2 for b in range(B - 1)])
        for k in lags:
            for b in range(W + k, B):
                Cp = np.einsum('bmn,bmk->nk', H[b - k - W:b - k].conj(), H[b - k - W:b - k])
                for r, store in ((2, pef2), (r95, pefe)):
                    P = topk_proj(Cp, r)
                    store[k][0] += np.linalg.norm(H[b] @ P) ** 2; store[k][1] += 1
            if k + W <= B:
                C0 = np.einsum('bmn,bmk->nk', H[:W].conj(), H[:W]); Ck = np.einsum('bmn,bmk->nk', H[k:k + W].conj(), H[k:k + W])
                V0, Vk = topk_proj(C0, r95), topk_proj(Ck, r95)
                ltov[k][0] += np.linalg.norm(V0.conj().T @ Vk) ** 2 / r95; ltov[k][1] += 1
    out.update({'r_eff90': acc['r90'] / n_runs, 'r_eff95': acc['r95'] / n_runs, 'participation_ratio': acc['pr'] / n_runs,
                'rho_S_inst_lag1': acc['rhoS1'] / n_runs,
                'PEF_r2': {str(k): pef2[k][0] / pef2[k][1] for k in lags},
                'PEF_reff95': {str(k): pefe[k][0] / pefe[k][1] for k in lags},
                'LT_overlap_reff95': {str(k): (ltov[k][0] / ltov[k][1] if ltov[k][1] else None) for k in lags}})
    print(f"{name:34s} D={Delta*1e3:5.1f}ms r90={out['r_eff90']:.1f} r95={out['r_eff95']:.1f} PR={out['participation_ratio']:.2f} "
          f"rhoS1={out['rho_S_inst_lag1']:.3f} PEF2(1)={out['PEF_r2']['1']:.3f} PEF2(16)={out['PEF_r2']['16']:.3f} "
          f"PEFe(1)={out['PEF_reff95']['1']:.3f} PEFe(48)={out['PEF_reff95']['48']:.3f} LT(16)={out['LT_overlap_reff95']['16']:.3f} LT(48)={out['LT_overlap_reff95']['48']:.3f}", flush=True)
    return out

if __name__ == '__main__':
    cases = {
      'NLOS_N32_R50_rho1.5_radial': dict(N=32, R=50.0, v=3.0, phi_v_deg=90.0),
      'NLOS_N32_R50_rho1.5_tangential': dict(N=32, R=50.0, v=3.0, phi_v_deg=0.0),
      'NLOS_N32_R50_rho0.5_tangential': dict(N=32, R=50.0, v=3.0, phi_v_deg=0.0, rho_loc=0.5),
      'LOS_K3dB_N32_R50_tangential': dict(N=32, R=50.0, v=3.0, phi_v_deg=0.0, K_los=2.0),
      'LOS_K10dB_N32_R50_tangential': dict(N=32, R=50.0, v=3.0, phi_v_deg=0.0, K_los=10.0),
      'NLOS_N16_R50_rho1.5_tangential': dict(N=16, R=50.0, v=3.0, phi_v_deg=0.0),
      'NLOS_N32_R100_rho1.5_tangential': dict(N=32, R=100.0, v=3.0, phi_v_deg=0.0),
    }
    res = {'seed': SEED, 'cases': {}}
    t0 = time.time()
    for name, kw in cases.items():
        for Delta in (2e-3, 10e-3):
            res['cases'][f"{name}__D{Delta*1e3:g}ms"] = run(name, kw, Delta)
            json.dump(res, open('t6_q12_longterm.json', 'w'), indent=1)
    print('done', time.time() - t0, 's')
