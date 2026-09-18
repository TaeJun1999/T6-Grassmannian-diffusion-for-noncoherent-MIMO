#!/usr/bin/env python3
"""
T6 gate-(0) ceiling pre-check (2026-09-16).

Model (assumptions A1-A6 in T6_ceiling_derivation.md):
    Y_b = X_b H_b + W_b,   X_b in C^{T x M},  H_b = G_b U_b^H in C^{M x N},
    U_b ~ Haar on V_r(C^N) (receive-side subspace S_b = span(H_b^H) = span(U_b)),
    G_b iid CN(0, gamma), W_b iid CN(0, 1), E||X_b||_F^2 = rho T  (rho = SNR).
Only M = r = 2 is supported (closed form F_D below is for rank-2 arguments).

Quantities (all per block, in nats; divide by T ln 2 for bit/s/Hz):
  * R_genie  : exact MI of the isotropic USTM input when U_b is known at the receiver
               (sufficient statistic Z = Y U, an iid noncoherent (T, M, r) channel at SNR gamma*rho).
  * GMI_svd  : generalized mutual information of the prior-free per-block SVD-projection receiver
               (U_hat = top-r right singular vectors of Y), a lower bound on the prior-free rate.
  * C_coh    : coherent Gaussian-input reference  T * E log det(I + rho/M G G^H).
  * finite K : exact MI (genie, structure-aware per-block ML via F_N, iid Rayleigh reference),
               GMI of GLRT-trace and SVD-projection metrics, and uncoded SER of every receiver.

Closed form (ours, verified by tests in test_t6_ceiling.py):
  F_D(b1, b2) = E_{U ~ Haar V_2(C^D)} exp(b1 ||U^H v1||^2 + b2 ||U^H v2||^2)
              = k^2 (k+1) [ e^{b1+b2} g_k(b1) g_k(b2) - (e^{b2} g_k(b2) - e^{b1} g_k(b1)) / (b2 - b1) ],
  k = D - 2,  g_k(b) = int_0^1 s^{k-1} e^{-b s} ds = gamma(k, b) / b^k.
"""
import argparse
import json
import time

import numpy as np
from scipy.special import gammainc, gammaln, logsumexp

LN2 = np.log(2.0)


# ----------------------------------------------------------------------------- sampling helpers
def crandn(rng, *shape):
    """iid CN(0,1)."""
    return (rng.standard_normal(shape) + 1j * rng.standard_normal(shape)) / np.sqrt(2.0)


def haar_stiefel(rng, n, D, r):
    """n Haar-distributed points on the complex Stiefel manifold V_r(C^D) (QR + phase fix, Mezzadri [M])."""
    A = crandn(rng, n, D, r)
    Q, R = np.linalg.qr(A)
    d = np.diagonal(R, axis1=-2, axis2=-1)
    ph = d / np.abs(d)
    return Q * ph[:, None, :]


def H_(A):
    """Batched conjugate transpose."""
    return np.conj(np.swapaxes(A, -1, -2))


# ----------------------------------------------------------------------------- special functions
def log_g(k, beta):
    """log g_k(beta),  g_k(beta) = int_0^1 s^{k-1} exp(-beta s) ds,  beta >= 0, k >= 1 (integer or real)."""
    beta = np.asarray(beta, dtype=float)
    out = np.empty_like(beta)
    small = beta < 1e-2
    if np.any(small):
        b = beta[small]
        acc = np.zeros_like(b)
        fact = 1.0
        for j in range(0, 14):
            if j > 0:
                fact *= j
            acc += (-b) ** j / (fact * (k + j))
        out[small] = np.log(acc)
    if np.any(~small):
        b = beta[~small]
        out[~small] = gammaln(k) + np.log(gammainc(k, b)) - k * np.log(b)
    return out


def log_h(k, beta):
    """log h_k(beta),  h_k(beta) = d/dbeta [e^beta g_k(beta)] = e^beta (g_k(beta) - g_{k+1}(beta)) > 0."""
    lgk = log_g(k, beta)
    lgk1 = log_g(k + 1, beta)
    return beta + lgk + np.log(-np.expm1(lgk1 - lgk))


_GL_NODES, _GL_WEIGHTS = np.polynomial.legendre.leggauss(6)     # on [-1, 1]
_GL_T = 0.5 * (_GL_NODES + 1.0)                                  # on [0, 1]
_GL_LOGW = np.log(0.5 * _GL_WEIGHTS)


def log_F(D, b1, b2, d_switch=0.5):
    """log F_D(b1, b2) for rank-2 arguments (b1, b2 >= 0), D >= 4. Vectorized over b1, b2 (broadcast).

    I2 = (e^{b2} g_k(b2) - e^{b1} g_k(b1)) / (b2 - b1) = int_0^1 h_k(b1 + t (b2 - b1)) dt.
    For |b1 - b2| <= d_switch the integral form is evaluated by 6-point Gauss-Legendre (no 1/d, no
    cancellation; relative error ~ d^12/(12!) ); otherwise the difference quotient is used (log domain).
    """
    b1 = np.asarray(b1, dtype=float)
    b2 = np.asarray(b2, dtype=float)
    b1, b2 = np.broadcast_arrays(b1, b2)
    k = D - 2
    if k < 2:
        raise ValueError("closed form needs D >= 4")
    logpref = 2.0 * np.log(k) + np.log(k + 1.0)
    with np.errstate(all="ignore"):
        L1 = b1 + log_g(k, b1)          # log( e^{b1} g_k(b1) )
        L2 = b2 + log_g(k, b2)
        logA = L1 + L2
        d = np.abs(b1 - b2)
        near = d <= d_switch
        # far branch: difference quotient
        Lmax = np.maximum(L1, L2)
        Lmin = np.minimum(L1, L2)
        logI2_far = Lmax + np.log(-np.expm1(-(Lmax - Lmin))) - np.log(d)
        # near branch: Gauss-Legendre on the integral representation
        bt = b1[..., None] + (b2 - b1)[..., None] * _GL_T          # (..., 6)
        logI2_near = logsumexp(_GL_LOGW + log_h(k, bt), axis=-1)
        logI2 = np.where(near, logI2_near, logI2_far)
        logF = logpref + logA + np.log(-np.expm1(logI2 - logA))
    return logF


def eig2x2_psd(p, q, z):
    """Eigenvalues (mu1 >= mu2 >= 0) of the Hermitian PSD matrix [[p, z], [conj(z), q]] (vectorized)."""
    half_tr = 0.5 * (p + q)
    disc = np.sqrt(np.maximum((0.5 * (p - q)) ** 2 + np.abs(z) ** 2, 0.0))
    mu1 = half_tr + disc
    det = np.maximum(p * q - np.abs(z) ** 2, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        mu2 = np.where(mu1 > 0, det / mu1, 0.0)
    return mu1, np.maximum(mu2, 0.0)


# ----------------------------------------------------------------------------- channel / signal
def snr_consts(rho, T, M, gamma):
    a = gamma * rho * T / M
    c = a / (1.0 + a)
    return a, c


def sample_channel(rng, n, T, M, N, r, rho, gamma, Phi):
    """Given input Phi (n, T, M) on V_M(C^T): returns X, G, U, Y = X G U^H + W and Z = Y U."""
    X = np.sqrt(rho * T / M) * Phi
    G = np.sqrt(gamma) * crandn(rng, n, M, r)
    U = haar_stiefel(rng, n, N, r)
    W = crandn(rng, n, T, N)
    Y = X @ G @ H_(U) + W
    Z = Y @ U
    return X, G, U, Y, Z


def svd_projection(Y, r):
    """Prior-free per-block estimate: U_hat = top-r right singular vectors of Y; returns Z_hat = Y U_hat and sv."""
    _, s, Vh = np.linalg.svd(Y, full_matrices=False)
    Uhat = H_(Vh[:, :r, :])
    return Y @ Uhat, s[:, :r]


# ----------------------------------------------------------------------------- continuous USTM input
def ustm_curves(T, M, N, r, snr_dB, n, gamma, rng, s_grid_u):
    assert M == 2 and r == 2, "closed form implemented for M = r = 2 only"
    rows = []
    for sdb in snr_dB:
        rho = 10 ** (sdb / 10.0)
        a, c = snr_consts(rho, T, M, gamma)
        Phi = haar_stiefel(rng, n, T, M)
        X, G, U, Y, Z = sample_channel(rng, n, T, M, N, r, rho, gamma, Phi)
        PhiH = H_(Phi)
        # --- genie (U known): exact per-block information density of the USTM input
        stat_g = np.sum(np.abs(PhiH @ Z) ** 2, axis=(1, 2))
        sv = np.linalg.svd(Z, compute_uv=False)
        info_g = c * stat_g - log_F(T, c * sv[:, 0] ** 2, c * sv[:, 1] ** 2)
        # --- exact upper bound on the gap:  Delta = I(X;U|Y) <= I(U;Y|X) = E[c||Phi^H Y U||^2 - log F_N(c mu1, c mu2)]
        A = PhiH @ Y                                             # (n, M, N)
        p = np.sum(np.abs(A[:, 0, :]) ** 2, axis=-1)
        q = np.sum(np.abs(A[:, 1, :]) ** 2, axis=-1)
        zz = np.sum(A[:, 0, :] * np.conj(A[:, 1, :]), axis=-1)
        mu1, mu2 = eig2x2_psd(p, q, zz)
        iuy = c * stat_g - log_F(N, c * mu1, c * mu2)
        # --- prior-free SVD-projection receiver: GMI over the scale s
        Zhat, sy = svd_projection(Y, r)
        stat_s = np.sum(np.abs(PhiH @ Zhat) ** 2, axis=(1, 2))
        s1, s2 = sy[:, 0] ** 2, sy[:, 1] ** 2
        best = (-np.inf, None, None)
        for u in s_grid_u:
            s = c * 10.0 ** u
            g = s * stat_s - log_F(T, s * s1, s * s2)
            gm = g.mean()
            if gm > best[0]:
                best = (gm, u, g.std(ddof=1) / np.sqrt(n))
        # --- coherent Gaussian-input reference
        GG = G @ H_(G)
        coh = T * np.linalg.slogdet(np.eye(M) + (rho / M) * GG)[1]
        rows.append(dict(
            snr_dB=float(sdb), a=a, c=c,
            R_genie=float(info_g.mean()), se_genie=float(info_g.std(ddof=1) / np.sqrt(n)),
            GMI_svd=float(best[0]), se_svd=float(best[2]), s_opt_u_svd=float(best[1]),
            C_coh=float(coh.mean()), se_coh=float(coh.std(ddof=1) / np.sqrt(n)),
            I_U_given_X=float(iuy.mean()), se_IUX=float(iuy.std(ddof=1) / np.sqrt(n)),
        ))
    return rows


# ----------------------------------------------------------------------------- finite constellations
def finiteK_curves(T, M, N, r, K, snr_dB, n, gamma, rng, s_grid_u, chunk=None):
    assert M == 2 and r == 2
    if chunk is None:
        chunk = int(np.clip(2 ** 23 // (K * N * 2), 1, 128))
    Phic = haar_stiefel(rng, K, T, M)          # random Haar Grassmannian constellation (fixed over SNR)
    PhicH = H_(Phic)                           # (K, M, T)
    S = len(s_grid_u)
    rows = []
    for sdb in snr_dB:
        rho = 10 ** (sdb / 10.0)
        a, c = snr_consts(rho, T, M, gamma)
        a0, c0 = snr_consts(rho, T, M, 1.0)    # iid Rayleigh reference (unit-variance entries)
        acc = {k_: [] for k_ in ["i_genie", "i_struct", "i_iid"]}
        gmi_glrt = np.zeros(S); gmi_svd = np.zeros(S)
        gmi_glrt_sq = np.zeros(S); gmi_svd_sq = np.zeros(S)
        err = dict(genie=0, struct=0, glrt=0, svd=0, iid=0)
        done = 0
        while done < n:
            b = min(chunk, n - done)
            jstar = rng.integers(0, K, b)
            Phi = Phic[jstar]
            X, G, U, Y, Z = sample_channel(rng, b, T, M, N, r, rho, gamma, Phi)
            H0 = crandn(rng, b, M, N)
            Y0 = X @ H0 + crandn(rng, b, T, N)
            # metrics, shape (b, K)
            A = np.einsum('kmt,btn->bkmn', PhicH, Y)
            p = np.sum(np.abs(A[:, :, 0, :]) ** 2, axis=-1)
            q = np.sum(np.abs(A[:, :, 1, :]) ** 2, axis=-1)
            zz = np.sum(A[:, :, 0, :] * np.conj(A[:, :, 1, :]), axis=-1)
            mu1, mu2 = eig2x2_psd(p, q, zz)
            l_glrt = mu1 + mu2
            l_struct = log_F(N, c * mu1, c * mu2)
            l_genie = c * np.sum(np.abs(np.einsum('kmt,btr->bkmr', PhicH, Z)) ** 2, axis=(2, 3))
            Zhat, _ = svd_projection(Y, r)
            l_svd = np.sum(np.abs(np.einsum('kmt,btr->bkmr', PhicH, Zhat)) ** 2, axis=(2, 3))
            l_iid = c0 * np.sum(np.abs(np.einsum('kmt,btn->bkmn', PhicH, Y0)) ** 2, axis=(2, 3))
            idx = np.arange(b)

            def dev(l):
                return l - l[idx, jstar][:, None]

            for name, l in [("genie", l_genie), ("struct", l_struct), ("iid", l_iid)]:
                acc["i_" + name].append(np.log(K) - logsumexp(dev(l), axis=1))
            for name, l, gm, gs in [("glrt", l_glrt, gmi_glrt, gmi_glrt_sq), ("svd", l_svd, gmi_svd, gmi_svd_sq)]:
                dv = dev(l)
                for si, u in enumerate(s_grid_u):
                    s = c * 10.0 ** u
                    g = np.log(K) - logsumexp(s * dv, axis=1)
                    gm[si] += g.sum(); gs[si] += (g ** 2).sum()
            for name, l in [("genie", l_genie), ("struct", l_struct), ("glrt", l_glrt), ("svd", l_svd), ("iid", l_iid)]:
                err[name] += int(np.sum(np.argmax(l, axis=1) != jstar))
            done += b
        row = dict(snr_dB=float(sdb), K=K, n=n, a=a, c=c)
        for name in ["genie", "struct", "iid"]:
            v = np.concatenate(acc["i_" + name])
            row["I_" + name] = float(v.mean()); row["se_" + name] = float(v.std(ddof=1) / np.sqrt(n))
        for name, gm, gs in [("glrt", gmi_glrt, gmi_glrt_sq), ("svd", gmi_svd, gmi_svd_sq)]:
            mean = gm / n
            si = int(np.argmax(mean))
            var = gs[si] / n - mean[si] ** 2
            row["GMI_" + name] = float(mean[si]); row["se_" + name] = float(np.sqrt(max(var, 0) / n))
            row["s_opt_u_" + name] = float(s_grid_u[si])
        for name in err:
            row["SER_" + name] = err[name] / n; row["nerr_" + name] = err[name]
        rows.append(row)
    return rows


def interp_cross(snr, ser, target=1e-2):
    """SNR (dB) at which log10(SER) crosses log10(target), piecewise-linear in dB; NaN if no crossing."""
    snr = np.asarray(snr, float); ser = np.asarray(ser, float)
    ls = np.log10(np.maximum(ser, 1e-12)); lt = np.log10(target)
    for i in range(len(snr) - 1):
        if (ls[i] - lt) * (ls[i + 1] - lt) <= 0 and ls[i] != ls[i + 1]:
            return float(snr[i] + (lt - ls[i]) * (snr[i + 1] - snr[i]) / (ls[i + 1] - ls[i]))
    return float("nan")


# ----------------------------------------------------------------------------- driver
def parse_snr(spec):
    a, b, c = [float(x) for x in spec.split(":")]
    return list(np.arange(a, b + 1e-9, c))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--configs", default="8,2,16,2;16,2,32,2")
    ap.add_argument("--snr", default="-10:20:2")
    ap.add_argument("--n_ustm", type=int, default=20000)
    ap.add_argument("--n_fin", type=int, default=4000)
    ap.add_argument("--Ks", default="64,1024")
    ap.add_argument("--gamma_mode", default="N/r", choices=["N/r", "1"])
    ap.add_argument("--K_low", type=int, default=0, help="extra large constellation evaluated on --snr_low only (0 = off)")
    ap.add_argument("--snr_low", default="-10:0:2")
    ap.add_argument("--n_low", type=int, default=600)
    ap.add_argument("--seed", type=int, default=20260916)
    ap.add_argument("--out", default="t6_ceiling_results.json")
    args = ap.parse_args()

    snr_dB = parse_snr(args.snr)
    s_grid_u = list(np.linspace(-3.0, 0.5, 36))
    results = dict(meta=dict(date="2026-09-16", seed=args.seed, snr_dB=snr_dB, n_ustm=args.n_ustm,
                             n_fin=args.n_fin, Ks=args.Ks, gamma_mode=args.gamma_mode,
                             numpy=np.__version__, s_grid_u=s_grid_u), configs=[])
    t0 = time.time()
    for ci, cfg in enumerate(args.configs.split(";")):
        T, M, N, r = [int(x) for x in cfg.split(",")]
        gamma = N / r if args.gamma_mode == "N/r" else 1.0
        rng = np.random.default_rng([args.seed, ci])
        entry = dict(T=T, M=M, N=N, r=r, gamma=gamma)
        print(f"[cfg] T={T} M={M} N={N} r={r} gamma={gamma}", flush=True)
        entry["ustm"] = ustm_curves(T, M, N, r, snr_dB, args.n_ustm, gamma, rng, s_grid_u)
        print(f"   ustm done  ({time.time()-t0:.0f}s)", flush=True)
        entry["finite"] = {}
        for K in [int(k) for k in args.Ks.split(",")]:
            nK = args.n_fin if K <= 256 else max(args.n_fin // 2, 200)
            entry["finite"][str(K)] = finiteK_curves(T, M, N, r, K, snr_dB, nK, gamma, rng, s_grid_u)
            rows = entry["finite"][str(K)]
            entry["finite"][str(K) + "_snr_at_ser1e-2"] = {
                name: interp_cross([x["snr_dB"] for x in rows], [x["SER_" + name] for x in rows])
                for name in ["genie", "struct", "glrt", "svd", "iid"]}
            print(f"   K={K} done  ({time.time()-t0:.0f}s)", flush=True)
        if args.K_low > 0:
            K = args.K_low
            entry["finite"][str(K)] = finiteK_curves(T, M, N, r, K, parse_snr(args.snr_low), args.n_low, gamma, rng, s_grid_u)
            print(f"   K_low={K} done  ({time.time()-t0:.0f}s)", flush=True)
        results["configs"].append(entry)
    with open(args.out, "w") as f:
        json.dump(results, f, indent=1)
    print("saved", args.out, f"({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
