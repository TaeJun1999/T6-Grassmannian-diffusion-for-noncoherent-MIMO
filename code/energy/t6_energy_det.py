#!/usr/bin/env python3
"""
t6_energy_det.py -- Energy-based (MCG 2016) detection in the T6 structured channel.

Model (M = 1, block length T, N receive antennas, rank-r receive-side subspace):
    Y = x h^T + W,  x = sqrt(p_k) 1_T,  h = U g,  U in V_r(C^N) (Haar),
    g ~ CN(0, (N/r) I_r)  (D-02: gamma = N/r, E||h||^2 = N),  W iid CN(0,1).
Energy receiver (MCG eq. (2), block-averaged):  S = ||Y||_F^2 / (N T).

Exact decomposition [정확]:
    N T S = c_k A + B,   c_k = p_k T N / r + 1,
    A ~ Gamma(r, 1),  B ~ Gamma(N T - r, 1),  A, B independent.
    E[S] = p_k + 1,   Var[S] = p_k^2 / r + 2 p_k /(N T) + 1/(N T).
iid (MCG) channel = r = N (gamma = 1).

Subspace-conditioned ML for OOK with known U (T6 genie, same input) [정확]:
    statistic ||U^H Y 1_T||^2 / T = c_k A   (A ~ Gamma(r,1)),  i.e. B removed.

Usage:
    python3 t6_energy_det.py --selftest
    python3 t6_energy_det.py --out t6_energy_det.json
"""
import argparse, json, sys
import numpy as np
from scipy import stats, integrate, optimize

SEED = 20260917


def _min_scalar(fun, lo, hi, ngrid=400):
    """Grid search + bounded Brent refinement (the SER curves are flat (=0.5) at both ends,
    which makes a plain golden-section search converge to a boundary)."""
    xs = np.linspace(lo, hi, ngrid)
    vals = np.array([fun(x) for x in xs])
    i = int(np.argmin(vals))
    a, b = xs[max(i - 1, 0)], xs[min(i + 1, ngrid - 1)]
    res = optimize.minimize_scalar(fun, bounds=(a, b), method="bounded", options={"xatol": (b - a) * 1e-4})
    if res.fun <= vals[i]:
        return float(res.fun), float(res.x)
    return float(vals[i]), float(xs[i])


def ser_binary_levels(p1, p2, N, T, r):
    """Exact SER of the energy receiver for two power levels p1 < p2 (equiprobable),
    with the SER-optimal threshold. p1 may be 0 (OOK). Returns (ser, tau)."""
    NT = N * T
    c1 = p1 * T * N / r + 1.0
    c2 = p2 * T * N / r + 1.0
    fA = stats.gamma(r)
    FB = stats.gamma(NT - r).cdf
    upper = fA.ppf(1 - 1e-14)

    def cdf_sum(tau, c):  # P(c A + B <= tau)
        if NT - r == 0:  # B degenerate (T = 1, r = N): S = c A / N
            return fA.cdf(tau / c)
        if p1 == 0 and c == c1:  # then c1 = 1 and A + B ~ Gamma(NT)
            return stats.gamma(NT).cdf(tau)
        hi = min(tau / c, upper)
        # split at the a where tau - c a crosses the bulk of B (mean +- 6 sd) to help quad
        mB, sB = NT - r, np.sqrt(NT - r)
        pts = sorted({min(max((tau - mB - k * sB) / c, 0.0), hi) for k in (-6, -3, 0, 3, 6)})
        pts = [p for p in pts if 0.0 < p < hi]
        val, _ = integrate.quad(lambda a: FB(tau - c * a) * fA.pdf(a), 0, hi, points=pts or None, limit=400)
        return val

    def ser(tau):
        return 0.5 * ((1 - cdf_sum(tau, c1)) + cdf_sum(tau, c2))

    mean1, mean2 = c1 * r + (NT - r), c2 * r + (NT - r)
    fun, tau = _min_scalar(ser, mean1 * 0.5, mean2, ngrid=60)
    return fun, tau / NT


def ser_floor_binary(kappa, r):
    """N -> infinity floor for two nonzero levels with ratio kappa = p2/p1 (r fixed):
    S -> p_k A / r + 1,  A ~ Gamma(r,1). Threshold on A: err = 0.5[P(A > t) + P(kappa^-1 ... )]."""
    F = stats.gamma(r).cdf

    def err(t):  # decide level 2 iff p A/r > t  <=>  A > t r / p ; use p1 = 1
        return 0.5 * ((1 - F(t)) + F(t / kappa))

    return _min_scalar(err, 1e-6, 30 * r)[0]


def ser_ook_subspace_ml(p, N, T, r):
    """OOK with known U: statistic c A vs A, A ~ Gamma(r,1), c = p T N / r + 1. Exact SER, optimal threshold."""
    c = p * T * N / r + 1.0
    F = stats.gamma(r).cdf

    def err(t):
        return 0.5 * ((1 - F(t)) + F(t / c))

    return _min_scalar(err, 1e-6, 30 * r * c)[0]


def snr_at_target(fun, target, lo_db, hi_db):
    """Find SNR (dB) at which fun(snr_dB) = target, fun decreasing in SNR. Returns None if not reached."""
    if fun(hi_db) > target:
        return None
    if fun(lo_db) < target:
        return None
    return float(optimize.brentq(lambda s: fun(s) - target, lo_db, hi_db, xtol=0.02))


def selftest():
    rng = np.random.default_rng(SEED)
    N, T, r, p = 32, 4, 2, 0.5
    n = 20000
    # Haar U
    Z = rng.standard_normal((N, r)) + 1j * rng.standard_normal((N, r))
    U, _ = np.linalg.qr(Z)
    S = np.empty(n)
    SU = np.empty(n)
    for i in range(n):
        g = np.sqrt(N / r / 2) * (rng.standard_normal(r) + 1j * rng.standard_normal(r))
        h = U @ g
        W = np.sqrt(0.5) * (rng.standard_normal((T, N)) + 1j * rng.standard_normal((T, N)))
        Y = np.sqrt(p) * np.ones((T, 1)) @ h[None, :] + W
        S[i] = np.linalg.norm(Y) ** 2 / (N * T)
        SU[i] = np.linalg.norm(U.conj().T @ Y.T @ np.ones(T)) ** 2 / T
    mean_th, var_th = p + 1, p ** 2 / r + 2 * p / (N * T) + 1 / (N * T)
    ok = True
    e1 = abs(S.mean() - mean_th) / mean_th
    e2 = abs(S.var() - var_th) / var_th
    print(f"[T1] E[S]  mc={S.mean():.4f} th={mean_th:.4f} relerr={e1:.3e}")
    print(f"[T2] Var[S] mc={S.var():.4f} th={var_th:.4f} relerr={e2:.3e}")
    ok &= e1 < 0.01 and e2 < 0.05
    c = p * T * N / r + 1
    m_th = c * r
    e3 = abs(SU.mean() - m_th) / m_th
    print(f"[T3] E[||U^H Y 1||^2/T] mc={SU.mean():.3f} th={m_th:.3f} relerr={e3:.3e}")
    ok &= e3 < 0.02
    # SER MC check for OOK energy receiver
    ser_th, tau = ser_binary_levels(0.0, 2 * p, N, T, r)
    nmc = 40000
    err = 0
    for i in range(nmc):
        on = i % 2 == 0
        g = np.sqrt(N / r / 2) * (rng.standard_normal(r) + 1j * rng.standard_normal(r))
        W = np.sqrt(0.5) * (rng.standard_normal((T, N)) + 1j * rng.standard_normal((T, N)))
        Y = (np.sqrt(2 * p) if on else 0.0) * np.ones((T, 1)) @ (U @ g)[None, :] + W
        s = np.linalg.norm(Y) ** 2 / (N * T)
        err += (s > tau) != on
    ser_mc = err / nmc
    print(f"[T4] OOK energy SER: exact={ser_th:.4f} mc={ser_mc:.4f} (tau={tau:.3f})")
    ok &= abs(ser_mc - ser_th) < 4 * np.sqrt(ser_th * (1 - ser_th) / nmc) + 1e-3
    # iid check: r = N reproduces MCG hardening variance (p+1)^2/N at T=1
    var_iid = p ** 2 / N + 2 * p / N + 1 / N
    print(f"[T5] iid variance formula: {var_iid:.5f} vs (p+1)^2/N={(p+1)**2/N:.5f}")
    ok &= abs(var_iid - (p + 1) ** 2 / N) < 1e-12
    print("selftest ok" if ok else "selftest FAILED")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--out", default="t6_energy_det.json")
    args = ap.parse_args()
    if args.selftest:
        sys.exit(0 if selftest() else 1)

    out = {"seed": SEED, "model": "M=1, x=sqrt(p_k) 1_T, h=U g, g~CN(0,(N/r)I), W~CN(0,1); S=||Y||_F^2/(NT)",
           "floor_binary": {}, "ook": {}, "ook_subspace_ml": {}}

    # 1) N -> infinity floors for two nonzero levels, ratio kappa, rank r
    for r in [1, 2, 4, 8, 32]:
        out["floor_binary"][str(r)] = {str(k): ser_floor_binary(k, r) for k in [2, 4, 8, 16, 64]}
    print("Floor SER (N->inf), two nonzero power levels with ratio kappa, rank r:")
    print("r \\ kappa  " + "  ".join(f"{k:>7d}" for k in [2, 4, 8, 16, 64]))
    for r, row in out["floor_binary"].items():
        print(f"r={r:>3s}     " + "  ".join(f"{v:7.4f}" for v in row.values()))

    # 2) OOK {0, 2 rho} energy receiver: SER vs SNR, SNR at 1e-2; structured (r=2) vs iid (r=N)
    target = 1e-2
    snr_grid = np.arange(-10, 21, 2.0)
    for (N, T, r) in [(32, 1, 2), (32, 16, 2), (32, 1, 32), (32, 16, 32), (16, 1, 2), (16, 8, 2), (16, 1, 16), (16, 8, 16),
                      (64, 1, 2), (64, 16, 2), (64, 1, 64), (64, 16, 64)]:
        key = f"N{N}_T{T}_r{r}"
        f = lambda s: ser_binary_levels(0.0, 2 * 10 ** (s / 10), N, T, r)[0]
        curve = {f"{s:.0f}": f(s) for s in snr_grid}
        snr1e2 = snr_at_target(f, target, -20, 25)
        out["ook"][key] = {"N": N, "T": T, "r": r, "rate_bpcu": 1 / T, "ser_vs_snr_dB": curve, "snr_dB_at_1e-2": snr1e2}
        print(f"OOK energy {key}: SNR@1e-2 = {snr1e2 if snr1e2 is None else round(snr1e2, 2)} dB")

    # 2b) prior-free block ML for OOK = matched-filter energy ||Y^T 1_T||^2 / T  (T=1 formula with p -> pT) [exact]
    out["ook_mf_priorfree_ml"] = {}
    for (N, T, r) in [(32, 16, 2), (16, 8, 2), (64, 16, 2)]:
        key = f"N{N}_T{T}_r{r}"
        f = lambda s: ser_binary_levels(0.0, 2 * 10 ** (s / 10) * T, N, 1, r)[0]
        snr1e2 = snr_at_target(f, target, -20, 25)
        out["ook_mf_priorfree_ml"][key] = {"N": N, "T": T, "r": r, "snr_dB_at_1e-2": snr1e2}
        print(f"OOK MF prior-free ML {key}: SNR@1e-2 = {snr1e2 if snr1e2 is None else round(snr1e2, 2)} dB")

    # 3) OOK with known subspace, ML (T6 genie, same input): SNR at 1e-2
    for (N, T, r) in [(32, 1, 2), (32, 16, 2), (16, 1, 2), (16, 8, 2), (64, 1, 2), (64, 16, 2)]:
        key = f"N{N}_T{T}_r{r}"
        f = lambda s: ser_ook_subspace_ml(2 * 10 ** (s / 10), N, T, r)
        snr1e2 = snr_at_target(f, target, -30, 25)
        out["ook_subspace_ml"][key] = {"N": N, "T": T, "r": r, "snr_dB_at_1e-2": snr1e2}
        print(f"OOK subspace-ML {key}: SNR@1e-2 = {snr1e2 if snr1e2 is None else round(snr1e2, 2)} dB")

    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
