#!/usr/bin/env python3
"""Unit tests for t6_ceiling.py (run: python3 test_t6_ceiling.py). All tests print PASS/FAIL and a summary."""
import numpy as np
from scipy.integrate import quad, dblquad
from scipy.special import logsumexp
import t6_ceiling as tc

rng = np.random.default_rng(1)
results = []


def check(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + ("  " + detail if detail else ""), flush=True)


# T1: F_D(beta, 0) equals the MGF of Beta(2, D-2) (numerical quadrature)
def t1():
    worst = 0.0
    for D in [4, 8, 16, 32]:
        k = D - 2
        for beta in [0.0, 0.3, 3.0, 30.0, 300.0]:
            f = lambda s: k * (k + 1) * s * (1 - s) ** (k - 1) * np.exp(beta * (s - 1))
            val, _ = quad(f, 0, 1, epsabs=0, epsrel=1e-12, limit=200)
            ref = beta + np.log(val)
            got = float(tc.log_F(D, beta, 0.0))
            worst = max(worst, abs(got - ref))
    check("T1 F_D(beta,0) vs Beta(2,D-2) MGF quadrature", worst < 1e-9, f"max|dlog|={worst:.2e}")


# T2: F_D(b1,b2) vs Haar Monte Carlo on V_2(C^D) (validates the truncated-Haar density claim [M])
def t2():
    worst_z = 0.0
    for D, b1, b2, n in [(8, 1.0, 0.5, 400000), (8, 3.0, 3.0, 400000), (8, 4.0, 1.0, 400000),
                         (16, 2.0, 1.0, 400000), (32, 2.5, 0.2, 400000), (4, 1.5, 0.7, 400000)]:
        U = tc.haar_stiefel(rng, n, D, 2)
        s1 = np.sum(np.abs(U[:, 0, :]) ** 2, axis=1)  # ||U^H e1||^2
        s2 = np.sum(np.abs(U[:, 1, :]) ** 2, axis=1)
        x = np.exp(b1 * s1 + b2 * s2)
        mc = x.mean(); se = x.std(ddof=1) / np.sqrt(n)
        cf = np.exp(float(tc.log_F(D, b1, b2)))
        z = abs(mc - cf) / se
        worst_z = max(worst_z, z)
    check("T2 F_D vs Haar MC (|z|<4)", worst_z < 4.0, f"max|z|={worst_z:.2f}")


# T2b: F_D vs direct 2-D quadrature of the derived marginal density (validates the algebra of the closed form)
def t2b():
    worst = 0.0
    for D, b1, b2 in [(8, 2.0, 5.0), (16, 10.0, 3.0), (8, 0.7, 0.2), (32, 20.0, 19.0)]:
        k = D - 2
        pref = k ** 2 * (k + 1)
        f = lambda s2, s1: pref * np.exp(b1 * s1 + b2 * s2 - b1 - b2) * (
            (1 - s1) ** (k - 1) * (1 - s2) ** (k - 1) - np.maximum(1 - s1 - s2, 0.0) ** (k - 1))
        val, _ = dblquad(f, 0, 1, 0, 1, epsabs=0, epsrel=1e-10)
        ref = b1 + b2 + np.log(val)
        got = float(tc.log_F(D, b1, b2))
        worst = max(worst, abs(got - ref))
    check("T2b F_D vs 2-D quadrature of marginal density", worst < 1e-7, f"max|dlog|={worst:.2e}")


# T3: full NT-dim Gaussian log-pdf (Jafar-Goldsmith eq.(12)-type covariance I + gamma conj(P) (x) X X^H)
#     equals the reduced log-pdf  -||Y||^2 + c||Phi^H Y U||^2 - TN log pi - M r log(1+a)
def t3():
    T, M, N, r = 8, 2, 16, 2
    worst = 0.0
    for rho in [0.05, 1.0, 30.0]:
        gamma = N / r
        a, c = tc.snr_consts(rho, T, M, gamma)
        Phi = tc.haar_stiefel(rng, 1, T, M)
        X, G, U, Y, Z = tc.sample_channel(rng, 1, T, M, N, r, rho, gamma, Phi)
        X0, U0, Y0 = X[0], U[0], Y[0]
        P = U0 @ U0.conj().T
        C = np.eye(N * T) + gamma * np.kron(np.conj(P), X0 @ X0.conj().T)  # vec column-stacking: n outer, t inner
        y = Y0.reshape(-1, order="F")
        full = -np.real(y.conj() @ np.linalg.solve(C, y)) - np.linalg.slogdet(C)[1] - N * T * np.log(np.pi)
        red = (-np.sum(np.abs(Y0) ** 2) + c * np.sum(np.abs(Phi[0].conj().T @ Y0 @ U0) ** 2)
               - T * N * np.log(np.pi) - M * r * np.log(1 + a))
        worst = max(worst, abs(full - red))
    check("T3 NT-dim JG-type pdf vs reduced pdf", worst < 1e-8, f"max|diff|={worst:.2e}")


# T4: structured prior-free marginal  E_U exp(c ||Phi^H Y U||^2)  (MC over U) vs  F_N(c mu1, c mu2)
def t4():
    T, M, N, r = 8, 2, 16, 2
    worst_z = 0.0
    for rho in [0.001, 0.01]:
        gamma = N / r
        a, c = tc.snr_consts(rho, T, M, gamma)
        Phi = tc.haar_stiefel(rng, 1, T, M)
        X, G, U, Y, Z = tc.sample_channel(rng, 1, T, M, N, r, rho, gamma, Phi)
        A = Phi[0].conj().T @ Y[0]                    # (M, N)
        B = A @ A.conj().T
        mu = np.linalg.eigvalsh(B)[::-1]
        cf = float(tc.log_F(N, c * mu[0], c * mu[1]))
        n = 400000
        Uu = tc.haar_stiefel(rng, n, N, r)
        x = np.exp(c * np.sum(np.abs(A[None] @ Uu) ** 2, axis=(1, 2)) - cf)
        z = abs(x.mean() - 1.0) / (x.std(ddof=1) / np.sqrt(n))
        worst_z = max(worst_z, z)
    check("T4 structured marginal closed form vs MC over U (|z|<4)", worst_z < 4.0, f"max|z|={worst_z:.2f}")


# T5: eig2x2_psd vs eigvalsh
def t5():
    A = tc.crandn(rng, 2000, 2, 5)
    B = A @ tc.H_(A)
    p, q, z = np.real(B[:, 0, 0]), np.real(B[:, 1, 1]), B[:, 0, 1]
    m1, m2 = tc.eig2x2_psd(p, q, z)
    ref = np.linalg.eigvalsh(B)
    err = max(np.max(np.abs(m1 - ref[:, 1])), np.max(np.abs(m2 - ref[:, 0])))
    check("T5 eig2x2 vs eigvalsh", err < 1e-9, f"max err={err:.2e}")


# T6: log_F vs 40-digit mpmath reference (both branches, incl. the switch at |b1-b2| = 0.5), symmetry, F(0,0)=1
def t6():
    import mpmath as mp
    mp.mp.dps = 40

    def g_mp(k, b):
        b = mp.mpf(b)
        return mp.mpf(1) / k if b == 0 else mp.gammainc(k, 0, b) / b ** k

    def logF_mp(D, b1, b2):
        k = D - 2; b1 = mp.mpf(b1); b2 = mp.mpf(b2)
        A = mp.e ** (b1 + b2) * g_mp(k, b1) * g_mp(k, b2)
        if b1 == b2:
            I2 = mp.e ** b1 * (g_mp(k, b1) - g_mp(k + 1, b1))
        else:
            I2 = (mp.e ** b2 * g_mp(k, b2) - mp.e ** b1 * g_mp(k, b1)) / (b2 - b1)
        return float(mp.log(k ** 2 * (k + 1) * (A - I2)))

    worst = 0.0; ok = True
    for D in [4, 8, 16, 32]:
        for b in [0.0, 0.005, 0.5, 5.0, 50.0, 5000.0]:
            for d in [0.0, 1e-9, 1e-6, 1e-3, 0.1, 0.49, 0.51, 2.0, 30.0]:
                ref = logF_mp(D, b, b + d)
                e = max(abs(float(tc.log_F(D, b, b + d)) - ref), abs(float(tc.log_F(D, b + d, b)) - ref))
                worst = max(worst, e)
        ok &= abs(float(tc.log_F(D, 0.0, 0.0))) < 1e-12
        seq = tc.log_F(D, np.linspace(0, 100, 50), 3.0)
        ok &= bool(np.all(np.diff(seq) > 0))
    ok &= worst < 1e-8
    check("T6 log_F vs mpmath reference / symmetry / F(0,0)=1 / monotone", ok, f"max|dlog|={worst:.2e}")


# T7: small finite-K run: exact-metric MI >= GMI(s) for all s; genie >= struct >= max(GMI_glrt, GMI_svd)
def t7():
    T, M, N, r, K = 8, 2, 16, 2, 64
    rows = tc.finiteK_curves(T, M, N, r, K, [0.0, 6.0], 1500, N / r, np.random.default_rng(7),
                             list(np.linspace(-3, 0.5, 36)), chunk=32)
    ok = True; det = []
    for row in rows:
        tol = 3 * (row["se_genie"] + row["se_struct"] + row["se_svd"] + row["se_glrt"])
        ok &= row["I_genie"] >= row["I_struct"] - tol
        ok &= row["I_struct"] >= max(row["GMI_svd"], row["GMI_glrt"]) - tol
        det.append((row["snr_dB"], row["I_genie"], row["I_struct"], row["GMI_svd"], row["GMI_glrt"]))
    # GMI(s) <= MI for the exact genie metric on the same data
    rng2 = np.random.default_rng(8)
    Phic = tc.haar_stiefel(rng2, K, T, M); n = 1500
    jstar = rng2.integers(0, K, n); rho = 1.0; gamma = N / r
    a, c = tc.snr_consts(rho, T, M, gamma)
    X, G, U, Y, Z = tc.sample_channel(rng2, n, T, M, N, r, rho, gamma, Phic[jstar])
    l = c * np.sum(np.abs(np.einsum('kmt,btr->bkmr', tc.H_(Phic), Z)) ** 2, axis=(2, 3))
    dv = l - l[np.arange(n), jstar][:, None]
    mi = np.log(K) - logsumexp(dv, axis=1)
    gm = [np.mean(np.log(K) - logsumexp(s * dv, axis=1)) for s in [0.1, 0.5, 1.0, 2.0, 10.0]]
    ok &= all(g <= mi.mean() + 1e-9 for g in gm) and abs(gm[2] - mi.mean()) < 1e-12
    check("T7 ordering genie>=struct>=GMIs and GMI(s)<=MI", ok, "; ".join(
        f"{d[0]:.0f}dB: {d[1]:.2f}/{d[2]:.2f}/{d[3]:.2f}/{d[4]:.2f}" for d in det))


if __name__ == "__main__":
    for t in [t1, t2, t2b, t3, t4, t5, t6, t7]:
        t()
    npass = sum(ok for _, ok in results)
    print(f"\n{npass}/{len(results)} tests passed")
