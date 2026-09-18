#!/usr/bin/env python3
"""Shared Q-12 metrics on a channel sequence H[b] in C^{M x N} (b = 0..B-1), used by the geometric model and the ray-tracing runs.
rho_G(k): |sum_b tr(G_b^H G~_{b+k})| / sqrt(...) with G_b = H_b U_b, G~_{b+k} = H_{b+k} U_b  (complex correlation magnitude, |alpha|)
rho_S(k): mean ||U_b^H U_{b+k}||_F^2 / r_inst  (instantaneous rank-r_inst subspace overlap)
r_eff90/95: eigenvalue count for 90/95 % energy of C = sum_b H_b^H H_b (whole sequence);  participation ratio
PEF_r(k):  E ||H_b P_past||^2 / ||H_b||^2, P_past = top-r eigenspace of the covariance of window [b-k-W, b-k)
LT_overlap(k): top-r_eff95 eigenspace overlap of windows [0,W) and [k,k+W)
"""
import numpy as np

def topk(C, r):
    ev, V = np.linalg.eigh(C); return V[:, -r:]

def metrics(H, r_inst=2, W=16, lags=(1, 2, 4, 8, 16, 32, 48)):
    H = np.asarray(H, complex); B, M, N = H.shape
    H = H / (np.linalg.norm(H, axis=(1, 2), keepdims=True) + 1e-30)
    lags = [k for k in lags if k < B]
    U = []
    for b in range(B):
        _, _, Vh = np.linalg.svd(H[b], full_matrices=False); U.append(Vh[:r_inst].conj().T)
    C = np.einsum('bmn,bmk->nk', H.conj(), H); ev = np.sort(np.linalg.eigvalsh(C))[::-1]; cs = np.cumsum(ev) / ev.sum()
    r90 = int(np.searchsorted(cs, 0.90) + 1); r95 = int(np.searchsorted(cs, 0.95) + 1)
    out = {'B': B, 'r_eff90': r90, 'r_eff95': r95, 'participation_ratio': float(ev.sum() ** 2 / (ev ** 2).sum()),
           'eig_top8_frac': [float(x) for x in ev[:8] / ev.sum()], 'rho_G': {}, 'rho_S': {}, 'PEF_r2': {}, 'PEF_reff95': {}, 'LT_overlap_reff95': {}}
    for k in lags:
        num = 0j; d1 = d2 = 0.0; ov = 0.0
        for b in range(B - k):
            G0 = H[b] @ U[b]; Gk = H[b + k] @ U[b]
            num += np.vdot(G0, Gk); d1 += np.linalg.norm(G0) ** 2; d2 += np.linalg.norm(Gk) ** 2
            ov += np.linalg.norm(U[b].conj().T @ U[b + k]) ** 2 / r_inst
        out['rho_G'][str(k)] = float(abs(num) / np.sqrt(d1 * d2)); out['rho_S'][str(k)] = float(ov / (B - k))
        p2 = pe = 0.0; n = 0
        for b in range(W + k, B):
            Cp = np.einsum('bmn,bmk->nk', H[b - k - W:b - k].conj(), H[b - k - W:b - k])
            p2 += np.linalg.norm(H[b] @ topk(Cp, 2)) ** 2; pe += np.linalg.norm(H[b] @ topk(Cp, r95)) ** 2; n += 1
        if n:
            out['PEF_r2'][str(k)] = float(p2 / n); out['PEF_reff95'][str(k)] = float(pe / n)
        if k + W <= B:
            C0 = np.einsum('bmn,bmk->nk', H[:W].conj(), H[:W]); Ck = np.einsum('bmn,bmk->nk', H[k:k + W].conj(), H[k:k + W])
            out['LT_overlap_reff95'][str(k)] = float(np.linalg.norm(topk(C0, r95).conj().T @ topk(Ck, r95)) ** 2 / r95)
    return out

if __name__ == '__main__':
    # selftest: (1) constant channel -> all ones; (2) iid Gaussian H -> rho_G ~ 0, rho_S ~ r/N, PEF_r2 ~ 2/N
    rng = np.random.default_rng(0); N, M, B = 32, 2, 80
    H0 = (rng.standard_normal((M, N)) + 1j * rng.standard_normal((M, N)))
    m = metrics(np.repeat(H0[None], B, 0)); print("[T1] constant:", m['rho_G']['16'], m['rho_S']['16'], m['PEF_r2']['16'], m['r_eff95'])
    ok = abs(m['rho_G']['16'] - 1) < 1e-9 and abs(m['rho_S']['16'] - 1) < 1e-9 and abs(m['PEF_r2']['16'] - 1) < 1e-9 and m['r_eff95'] == 2
    Hi = rng.standard_normal((B, M, N)) + 1j * rng.standard_normal((B, M, N))
    m = metrics(Hi); print("[T2] iid:", round(m['rho_G']['1'], 3), round(m['rho_S']['1'], 3), 'expect ~', round(2 / N, 3), round(m['PEF_r2']['1'], 3), m['r_eff95'])
    ok &= m['rho_G']['1'] < 0.15 and abs(m['rho_S']['1'] - 2 / N) < 0.04 and abs(m['PEF_r2']['1'] - 2 / N) < 0.04
    print("selftest ok" if ok else "selftest FAILED")
