#!/usr/bin/env python3
"""
T6 Phase 0 item 3 (Q-06): pre-estimate of what differential USTM gains on the structured channel.

Model (same as t6_ceiling.py, plus block-to-block correlation):
    Y_b = X_b H_b + W_b,   H_b = G_b U^H,   U ~ Haar V_r(C^N) FIXED over a run (subspace persists),
    G_b: AR(1) across blocks,  G_b = alpha G_{b-1} + sqrt(1 - alpha^2) E_b,  E_b iid CN(0, gamma),
    so E[G_b G_{b-k}^H] ~ alpha^k (alpha = 1: channel constant; alpha = 0: gains iid, only the subspace persists).
    Jakes correspondence (one block apart): alpha = J0(2 pi f_D T T_s)  [M, Jakes].

Differential USTM, tall-block generalization of Hochwald-Sweldens 2000 [V, web copy] (HS use T = M):
    constellation {V_z} of K Haar-random T x T unitaries;  X_0 = sqrt(rho T / M) Phi_0,  X_b = V_{z_b} X_{b-1};
    receiver (HS eq. (21)/(22)):  z_hat = argmin_z ||Y_b - V_z Y_{b-1}||_F = argmax_z Re tr(V_z^H Y_b Y_{b-1}^H).
    Variants: 'diff'   : raw Y (prior-free, N noise dimensions),
              'diffU'  : Y projected on the true U (genie subspace + differential),
              'diffS2' : Y projected on the top-r right singular subspace of [Y_{b-1}; Y_b] (prior-free 2-block SVD).
Per-block references on the same channel statistics: structure-aware ML ('struct', isotropic-prior marginalization
via F_N) and genie-U ML ('genie'), random-Haar K-point Grassmannian constellation, as in t6_ceiling.py.
Rate is matched: log2 K bits per block of T channel uses for every scheme; E||X_b||_F^2 = rho T.

Output: SER vs SNR for every scheme and the SNR at SER 1e-2 (interp_cross), per (config, K, alpha).
"""
import argparse, json, os, sys, time
import numpy as np
from scipy.special import j0, logsumexp
from scipy.optimize import brentq


def _import_ceiling(path):
    sys.path.insert(0, path)
    import t6_ceiling as C
    return C


def alpha_to_fdTs(alpha, T):
    """f_D T_s such that J0(2 pi f_D T T_s) = alpha (smallest root); alpha = 1 -> 0."""
    if alpha >= 1.0:
        return 0.0
    x = brentq(lambda x: j0(x) - alpha, 1e-9, 2.41)
    return x / (2 * np.pi * T)


def ar1_gains(rng, C, R, B, M, r, gamma, alpha):
    G = np.empty((R, B, M, r), dtype=complex)
    G[:, 0] = np.sqrt(gamma) * C.crandn(rng, R, M, r)
    s = np.sqrt(gamma * (1.0 - alpha ** 2))
    for b in range(1, B):
        G[:, b] = alpha * G[:, b - 1] + s * C.crandn(rng, R, M, r)
    return G


def diff_metric(Ycur, Yprev, VconjFlat, chunk=1000):
    """m[d, z] = Re tr(V_z^H Ycur_d Yprev_d^H); Ycur, Yprev: (D, T, n); VconjFlat: (K, T*T)."""
    D, T, _ = Ycur.shape
    out = np.empty((D, VconjFlat.shape[0]))
    for i in range(0, D, chunk):
        A = Ycur[i:i + chunk] @ np.conj(np.swapaxes(Yprev[i:i + chunk], -1, -2))     # (d, T, T)
        out[i:i + chunk] = np.real(A.reshape(A.shape[0], T * T) @ VconjFlat.T)
    return out


def sim_differential(C, T, M, N, r, K, rho, gamma, alpha, R, B, rng, V, VconjFlat, variants):
    """One SNR point. Returns dict variant -> (n_err, n_dec)."""
    G = ar1_gains(rng, C, R, B, M, r, gamma, alpha)          # (R, B, M, r)
    U = C.haar_stiefel(rng, R, N, r)                          # fixed per run
    UH = C.H_(U)
    z = rng.integers(0, K, (R, B))
    X = np.sqrt(rho * T / M) * C.haar_stiefel(rng, R, T, M)   # X_0
    Y = np.empty((R, B, T, N), dtype=complex)
    for b in range(B):
        if b > 0:
            X = V[z[:, b]] @ X
        Y[:, b] = X @ (G[:, b] @ UH) + C.crandn(rng, R, T, N)
    Yprev = Y[:, :-1].reshape(-1, T, N)
    Ycur = Y[:, 1:].reshape(-1, T, N)
    ztrue = z[:, 1:].reshape(-1)
    out = {}
    if 'diff' in variants:
        m = diff_metric(Ycur, Yprev, VconjFlat)
        out['diff'] = (int(np.sum(np.argmax(m, axis=1) != ztrue)), ztrue.size)
    if 'diffU' in variants:
        Urep = np.repeat(U[:, None], B - 1, axis=1).reshape(-1, N, r)
        m = diff_metric(Ycur @ Urep, Yprev @ Urep, VconjFlat)
        out['diffU'] = (int(np.sum(np.argmax(m, axis=1) != ztrue)), ztrue.size)
    if 'diffS2' in variants:
        S = np.concatenate([Yprev, Ycur], axis=1)              # (D, 2T, N)
        _, _, Vh = np.linalg.svd(S, full_matrices=False)
        Uhat = C.H_(Vh[:, :r, :])                              # (D, N, r)
        m = diff_metric(Ycur @ Uhat, Yprev @ Uhat, VconjFlat)
        out['diffS2'] = (int(np.sum(np.argmax(m, axis=1) != ztrue)), ztrue.size)
    return out


def sim_perblock(C, T, M, N, r, K, rho, gamma, n, rng, Phic, PhicH, chunk):
    """Per-block struct-ML and genie SER on independent blocks (their SER does not depend on alpha)."""
    a, c = C.snr_consts(rho, T, M, gamma)
    err = dict(struct=0, genie=0)
    done = 0
    while done < n:
        b = min(chunk, n - done)
        jstar = rng.integers(0, K, b)
        X, G, U, Y, Z = C.sample_channel(rng, b, T, M, N, r, rho, gamma, Phic[jstar])
        A = np.einsum('kmt,btn->bkmn', PhicH, Y)
        p = np.sum(np.abs(A[:, :, 0, :]) ** 2, axis=-1)
        q = np.sum(np.abs(A[:, :, 1, :]) ** 2, axis=-1)
        zz = np.sum(A[:, :, 0, :] * np.conj(A[:, :, 1, :]), axis=-1)
        mu1, mu2 = C.eig2x2_psd(p, q, zz)
        l_struct = C.log_F(N, c * mu1, c * mu2)
        l_genie = c * np.sum(np.abs(np.einsum('kmt,btr->bkmr', PhicH, Z)) ** 2, axis=(2, 3))
        err['struct'] += int(np.sum(np.argmax(l_struct, axis=1) != jstar))
        err['genie'] += int(np.sum(np.argmax(l_genie, axis=1) != jstar))
        done += b
    return {k: (v, n) for k, v in err.items()}


def selftest(C):
    rng = np.random.default_rng(1)
    T, N, K = 6, 5, 7
    V = C.haar_stiefel(rng, K, T, T)
    Vc = np.conj(V).reshape(K, T * T)
    Yp, Yc = C.crandn(rng, 3, T, N), C.crandn(rng, 3, T, N)
    m = diff_metric(Yc, Yp, Vc)
    ref = np.array([[-0.5 * (np.linalg.norm(Yc[d] - V[k] @ Yp[d]) ** 2 - np.linalg.norm(Yc[d]) ** 2 - np.linalg.norm(Yp[d]) ** 2)
                     for k in range(K)] for d in range(3)])
    assert np.allclose(m, ref, atol=1e-9), "HS identity failed"
    G = ar1_gains(rng, C, 4000, 3, 2, 2, 8.0, 0.9)
    var = np.mean(np.abs(G) ** 2)
    corr = np.mean(np.real(G[:, 1] * np.conj(G[:, 0]))) / var
    assert abs(var - 8.0) < 0.4 and abs(corr - 0.9) < 0.03, (var, corr)
    assert abs(alpha_to_fdTs(1.0, 16)) == 0.0 and abs(j0(2 * np.pi * 16 * alpha_to_fdTs(0.9, 16)) - 0.9) < 1e-9
    # alpha = 1, high SNR -> no errors; alpha = 0 -> chance
    V = C.haar_stiefel(rng, 16, 8, 8); Vc = np.conj(V).reshape(16, 64)
    o1 = sim_differential(C, 8, 2, 16, 2, 16, 10 ** 2.0, 8.0, 1.0, 100, 6, rng, V, Vc, ['diff'])
    o0 = sim_differential(C, 8, 2, 16, 2, 16, 10 ** 2.0, 8.0, 0.0, 100, 6, rng, V, Vc, ['diff'])
    assert o1['diff'][0] == 0, o1
    # alpha = 0 with a persistent subspace is NOT chance (Y_b Y_{b-1}^H ~ V_z X (G_b G_{b-1}^H) X^H keeps the
    # correct rotation's correlation large, with random sign): expect SER well above 0.5 but below chance.
    ser0 = o0['diff'][0] / o0['diff'][1]
    assert 0.5 < ser0 < 1 - 1 / 16, o0
    # alpha = 0 AND an iid subspace per block -> chance
    rng2 = np.random.default_rng(3)
    T_, N_, K_ = 8, 16, 16
    Yp, Yc = C.crandn(rng2, 400, T_, N_), C.crandn(rng2, 400, T_, N_)
    m = diff_metric(Yc, Yp, Vc)
    zt = rng2.integers(0, K_, 400)
    serc = np.mean(np.argmax(m, axis=1) != zt)
    assert abs(serc - (1 - 1 / K_)) < 0.06, serc
    print("selftest ok")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ceiling_dir', default='/mnt/project')
    ap.add_argument('--configs', default='8,2,16,2;16,2,32,2')
    ap.add_argument('--snr', default='-14:8:1', help='dB grid a:b:step (use --snr=-14:8:1)')
    ap.add_argument('--Ks', default='64,1024')
    ap.add_argument('--alphas', default='1.0,0.984,0.9,0.5,0.0')
    ap.add_argument('--alphas_proj', default='1.0,0.984,0.9', help='alphas for which diffU/diffS2 are also run')
    ap.add_argument('--R', type=int, default=200, help='runs per SNR point')
    ap.add_argument('--B', type=int, default=21, help='blocks per run (B-1 differential decisions)')
    ap.add_argument('--n_pb', type=int, default=4000, help='independent blocks for per-block references')
    ap.add_argument('--seed', type=int, default=20260916)
    ap.add_argument('--out', default='t6_diff_results.json')
    ap.add_argument('--selftest', action='store_true')
    args = ap.parse_args()
    C = _import_ceiling(args.ceiling_dir)
    if args.selftest:
        selftest(C); return
    snr_specs = args.snr.split(';')                      # one spec, or one per config
    alphas = [float(x) for x in args.alphas.split(',')]
    alphas_proj = [float(x) for x in args.alphas_proj.split(',')]
    res = dict(meta=dict(date='2026-09-16', seed=args.seed, snr=args.snr, R=args.R, B=args.B, n_pb=args.n_pb,
                         alphas=alphas, alphas_proj=alphas_proj, numpy=np.__version__), configs=[])
    t0 = time.time()
    for ci, cfg in enumerate(args.configs.split(';')):
        T, M, N, r = [int(x) for x in cfg.split(',')]
        gamma = N / r
        snr_dB = C.parse_snr(snr_specs[min(ci, len(snr_specs) - 1)])
        entry = dict(T=T, M=M, N=N, r=r, gamma=gamma, snr_dB=snr_dB, K={})
        print(f"[cfg] T={T} M={M} N={N} r={r} gamma={gamma}  alphas->fdTs: "
              + ", ".join(f"{al}:{alpha_to_fdTs(al, T):.4f}" for al in alphas), flush=True)
        for K in [int(k) for k in args.Ks.split(',')]:
            rng = np.random.default_rng([args.seed, ci, K])
            V = C.haar_stiefel(rng, K, T, T); VconjFlat = np.conj(V).reshape(K, T * T)
            Phic = C.haar_stiefel(rng, K, T, M); PhicH = C.H_(Phic)
            chunk = int(np.clip(2 ** 23 // (K * N * 2), 1, 128))
            n_pb = args.n_pb if K <= 256 else max(args.n_pb // 2, 200)   # same rule as t6_ceiling.py
            rows = []
            for sdb in snr_dB:
                rho = 10 ** (sdb / 10.0)
                row = dict(snr_dB=float(sdb))
                pb = sim_perblock(C, T, M, N, r, K, rho, gamma, n_pb, rng, Phic, PhicH, chunk)
                for k, (e, n) in pb.items():
                    row[f'SER_{k}'] = e / n; row[f'nerr_{k}'] = e; row[f'n_{k}'] = n
                for al in alphas:
                    variants = ['diff'] + (['diffU', 'diffS2'] if al in alphas_proj else [])
                    o = sim_differential(C, T, M, N, r, K, rho, gamma, al, args.R, args.B, rng, V, VconjFlat, variants)
                    for k, (e, n) in o.items():
                        row[f'SER_{k}_a{al}'] = e / n; row[f'nerr_{k}_a{al}'] = e; row[f'n_{k}_a{al}'] = n
                rows.append(row)
                print(f"   K={K} SNR={sdb:5.1f}  struct={row['SER_struct']:.4f} genie={row['SER_genie']:.4f}  "
                      + " ".join(f"diff(a{al})={row[f'SER_diff_a{al}']:.4f}" for al in alphas)
                      + f"  ({time.time()-t0:.0f}s)", flush=True)
            names = ['struct', 'genie'] + [f'diff_a{al}' for al in alphas] \
                + [f'{v}_a{al}' for al in alphas_proj for v in ('diffU', 'diffS2')]
            cross = {nm: C.interp_cross([x['snr_dB'] for x in rows], [x[f'SER_{nm}'] for x in rows]) for nm in names}
            entry['K'][str(K)] = dict(rows=rows, snr_at_ser1em2=cross)
            print(f"   K={K} crossings(dB): " + ", ".join(f"{k}={v:.2f}" for k, v in cross.items()), flush=True)
        res['configs'].append(entry)
        with open(args.out, 'w') as f:
            json.dump(res, f, indent=1)
    print("saved", args.out, f"({time.time()-t0:.0f}s)")


if __name__ == '__main__':
    main()
