#!/usr/bin/env python3
"""
Alpha-dependence of the strongest full-channel competitor (Q-11 (b), upper bound):
  predicted-CSI coherent detector = coherent ML with H_hat_b = alpha * H_{b-1}, where H_{b-1} is known EXACTLY (genie previous
  channel; real trackers are worse). Dynamics as in T6_diff_ustm_note.md §3: H_b = G_b U^H, U fixed Haar on V_r(C^N),
  G_b = alpha G_{b-1} + sqrt(1-alpha^2) E_b, E_b iid CN(0, N/r). Input: random-Haar Grassmannian constellation, K points,
  X_k = sqrt(rho T / M) Phi_k. Metric: argmin_k ||Y - X_k H_hat||_F^2  (mismatched coherent ML).
Reference points from t6_ceiling_results.json / t6_ceiling_LT.json: prior-free struct-ML, LT genie (r_eff), exact-U genie.
Output: SER vs SNR for alpha in list; SNR @ SER 1e-2. seed 20260917.
"""
import json, argparse, numpy as np, sys
sys.path.insert(0, '/mnt/project')
from t6_ceiling import haar_stiefel, crandn, interp_cross

def run(T, M, N, r, K, alphas, snrs, n, rng):
    gamma = N / r
    Phi = haar_stiefel(rng, K, T, M)                                   # (K, T, M)
    out = {}
    for alpha in alphas:
        sers = []
        for snr in snrs:
            rho = 10 ** (snr / 10); X = np.sqrt(rho * T / M) * Phi   # (K, T, M)
            err = 0
            for i in range(n):
                U = haar_stiefel(rng, 1, N, r)[0]                       # (N, r)
                Gp = np.sqrt(gamma) * crandn(rng, M, r)                 # previous block gain
                G = alpha * Gp + np.sqrt(1 - alpha ** 2) * np.sqrt(gamma) * crandn(rng, M, r)
                H = G @ U.conj().T; Hp = Gp @ U.conj().T; Hhat = alpha * Hp
                k = rng.integers(K)
                Y = X[k] @ H + crandn(rng, T, N)
                d = np.linalg.norm(Y[None] - X @ Hhat[None], axis=(1, 2)) ** 2   # (K,)
                err += int(np.argmin(d) != k)
            sers.append(err / n)
        sers = np.array(sers)
        snr_at = interp_cross(np.array(snrs, float), sers, 1e-2) if hasattr(sys.modules['t6_ceiling'], 'interp_cross') else None
        out[str(alpha)] = {'snr_dB': list(map(float, snrs)), 'ser': sers.tolist(), 'snr_at_ser1e-2': None if snr_at is None else float(snr_at)}
        print(f"alpha={alpha}: SER {np.round(sers,4)}  SNR@1e-2 = {out[str(alpha)]['snr_at_ser1e-2']}", flush=True)
    return out

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--n', type=int, default=3000); ap.add_argument('--out', default='t6_alpha_competitor.json')
    a = ap.parse_args(); rng = np.random.default_rng(20260917)
    res = {'seed': 20260917, 'model': 'predicted-CSI coherent ML, H_hat = alpha H_{b-1} (H_{b-1} exact)', 'cases': {}}
    for (T, M, N, r, K, snrs) in [(16, 2, 32, 2, 64, list(range(-14, 1, 1))), (8, 2, 16, 2, 64, list(range(-8, 7, 1)))]:
        print(f"== ({T},{M},{N},{r}) K={K}", flush=True)
        res['cases'][f"{T},{M},{N},{r},K{K}"] = run(T, M, N, r, K, [1.0, 0.95, 0.9, 0.8, 0.6, 0.0], snrs, a.n, rng)
        json.dump(res, open(a.out, 'w'), indent=1)
    print('done')
