#!/usr/bin/env python3
"""
T6 Q-15 (2026-09-18): detector comparison on the ray-tracing channel sequences.

First direct evidence for gate (i): how much room is left between the best realistic
baseline and the long-term-subspace ceiling, on the REAL RT dynamics (not a synthetic AR(1)).

Model per block b (handoff notation):
    Y_b = X_b H_b + W_b,   X_b = sqrt(rho T / M) Phi_{k_b},  Phi_k^H Phi_k = I_M,
    H_b in C^{M x N} from ray tracing (t6_q12_rt_H*.json),  W_b iid CN(0, 1).
Normalization (matches t6_ceiling.py: E||H||_F^2 = M N, rho = per-rx-antenna average SNR):
    --norm sequence  : one scale per trajectory, mean_b ||H_b||_F^2 = M N   (DEFAULT; keeps
                       the block-to-block power fluctuation, i.e. the small-scale fading a
                       predictor can exploit and a noncoherent receiver suffers from)
    --norm block     : per-block scale (removes amplitude dynamics; sensitivity check only)

Unified likelihood (derivation in the session notes, EXACT under its own model):
  If the receiver models H = G Lam^{1/2} Q^H with G iid CN(0,1), sum_i lam_i = N, then the
  columns of Y Q are independent and
      log p(Y | k) = const + sum_i c_i || Phi_k^H Y q_i ||^2,  c_i = a_i/(1+a_i), a_i = lam_i rho T/M.
  lam = 1 (all i)                    -> prior-free GLRT
  lam = N/r on span(U_b), 0 else     -> exact-U genie
  lam = N/r_eff on U_LT, 0 else      -> long-term-subspace projected conditional ML
  lam = estimated eigenspectrum      -> tracker + soft conditional ML (>= hard projection)
Outside that family:
  struct-ML : rank-r isotropic-prior marginalization, log F_N(c mu1, c mu2)  (t6_ceiling.log_F)
  error-aware coherent : H = Hhat + E, E iid CN(0, sig_e^2) ->
      metric_k = -||Y - X_k Hhat||_F^2 + c_e ||Phi_k^H (Y - X_k Hhat)||_F^2
  2T-block USTM : Psi_{k1k2} = [Phi_k1; Phi_k2]/sqrt(2) in V_M(C^{2T}), prior-free struct-ML at T->2T.
All of these are EXACT for their assumed model and APPROXIMATE on RT channels (mismatch).

Detector names (--detectors, comma separated; 'all' = DEFAULT_SET):
  glrt      prior-free GLRT                                   (no side info)
  struct    prior-free structure-aware ML (D-12 baseline #1)   (no side info)
  ltO<r>    oracle-window LT tracker, hard projection + cond. ML, rank r     [upper bound on causal tracking]
  ltOs<r>   same but rank-2 isotropic-in-subspace struct-ML (needs r >= 4)
  ltD<r>    decision-directed LT tracker, hard projection + cond. ML, rank r [realistic]
  covO      oracle-window full-spectrum soft conditional ML
  covD      decision-directed full-spectrum soft conditional ML              [realistic, usually strongest tracker]
  predO     genie-past AR(1) full-channel predictor + error-aware coherent   [upper bound on Q-11(b)]
  predD     decision-directed AR(1) predictor + error-aware coherent         [realistic]
  ustm2T    2T-block USTM (channel constant over 2 blocks), K^2 points
  genieU    exact instantaneous U_b (gate-(0) exact-U genie)
  genieLT<r>  exact NON-CAUSAL long-term subspace of the WHOLE trajectory, rank r (static LT genie)
  genieLTw<r> exact NON-CAUSAL long-term subspace of a CENTRED window of W blocks around b
              -> the proper ceiling for a perfect prior that TRACKS S_b^LT (it drifts: LT(32)=0.71-0.75
              on the RT NLOS trajectories), and therefore the ceiling gate (i) must be read against
  genieCovw   centred-window full-spectrum soft ML: full second-order knowledge (report only; this is
              more than the D-14 prior object, which is the subspace)
  genieH    coherent genie (exact H_b)

Output: per detector the SER curve, the SNR at SER 1e-2, a marginal cluster-bootstrap CI, and the
per-block error counts err_block (S x B) -- the last one is what t6_q15_summarize.py needs for the
PAIRED cluster bootstrap on detector DIFFERENCES (the gate reads differences, not absolute levels).

Statistics: B blocks x L noise repetitions per SNR. Errors within a block share H_b, so the
confidence interval is a CLUSTER bootstrap over blocks, not a binomial CI. SNR at SER 1e-2 by
t6_ceiling.interp_cross. All detectors see the SAME symbols and the SAME noise (paired comparison
-> the dB GAPS are far more accurate than the individual curves).

Cost (1 core, B=80, L=200, 8 SNR points, one trajectory): ~10 min for everything except ustm2T,
which costs ~10x more (K^2 = 4096 log_F evaluations per pair) -> run it in its own job.

Usage:
    python3 t6_q15_detectors.py --selftest                      # 2 validations, ~5 min
    python3 t6_q15_detectors.py --H t6_q12_rt_H_q15.json --detectors glrt,struct,ltO4,ltO8,ltD4,ltD8,covO,covD,predO,predD,genieU,genieLT4,genieLT8,genieH --out t6_q15_main.json
    python3 t6_q15_detectors.py --H t6_q12_rt_H_q15.json --detectors ustm2T --out t6_q15_ustm2T.json
    python3 t6_q15_detectors.py --H ... --quick                 # smoke test
Requires t6_ceiling.py on the path (log_F, eig2x2_psd, haar_stiefel, crandn, interp_cross).
"""
import argparse, base64, json, os, re, sys, time, zlib
import numpy as np

sys.path[:0] = [os.path.dirname(os.path.abspath(__file__)), os.getcwd(),
                os.path.dirname(os.path.dirname(os.path.abspath(__file__)))]
import t6_ceiling as C      # must be reachable: keep it next to this file or set PYTHONPATH

DEFAULT_SET = ("glrt,struct,ltO4,ltO8,ltD4,ltD8,covO,covD,predO,predD,ustm2T,"
               "genieU,genieLT4,genieLT8,genieLTw4,genieLTw8,genieCovw,genieH")


# --------------------------------------------------------------------------- data
def load_H(paths):
    runs = {}
    for p in paths:
        d = json.load(open(p))
        for k, v in d["runs"].items():
            H = np.frombuffer(base64.b64decode(v["b64"]), dtype=v["dtype"]).reshape(v["shape"])
            runs[k] = H.astype(np.complex128)
    return runs


def normalize(H, mode="sequence"):
    M, N = H.shape[1], H.shape[2]
    if mode == "sequence":
        s = np.mean(np.linalg.norm(H, axis=(1, 2)) ** 2)
        return H * np.sqrt(M * N / s)
    if mode == "block":
        s = np.linalg.norm(H, axis=(1, 2), keepdims=True) ** 2
        return H * np.sqrt(M * N / s)
    raise ValueError(mode)


def topk_eig(Cmat, r):
    """top-r eigenvectors (columns) of a Hermitian matrix; batched over leading axes."""
    _, V = np.linalg.eigh(Cmat)
    return V[..., -r:]


# --------------------------------------------------------------------------- metrics
def eig2(A):
    """mu1 >= mu2: eigenvalues of A A^H for A with M = 2 rows. A: (..., 2, D)."""
    p = np.sum(np.abs(A[..., 0, :]) ** 2, axis=-1)
    q = np.sum(np.abs(A[..., 1, :]) ** 2, axis=-1)
    z = np.sum(A[..., 0, :] * np.conj(A[..., 1, :]), axis=-1)
    return C.eig2x2_psd(p, q, z)


def A_of(PhiH, Y):
    """Phi_k^H Y for every codeword. PhiH: (K, M, T), Y: (L, T, D) -> (L, K, M, D)."""
    return np.einsum("kmt,ltd->lkmd", PhiH, Y)


def m_glrt(A):
    mu1, mu2 = eig2(A)
    return mu1 + mu2


def m_struct(A, D, c):
    mu1, mu2 = eig2(A)
    return C.log_F(D, c * mu1, c * mu2)


def m_cov(A, cvec):
    """sum_i c_i ||Phi_k^H Y q_i||^2 ; A = Phi^H (Y Q) : (L, K, M, D), cvec: (D,) or (L, D)."""
    e = np.sum(np.abs(A) ** 2, axis=2)                      # (L, K, D)
    if cvec.ndim == 1:
        return e @ cvec
    return np.einsum("lkd,ld->lk", e, cvec)


def m_coh(Y, PhiH, Hhat, rho, T, M, sig2e):
    """Error-aware coherent: -||Y - X_k Hhat||^2 + c_e ||Phi_k^H (Y - X_k Hhat)||^2.
    Y: (L, T, N), Hhat: (L, M, N) or (M, N), sig2e: (L,) or scalar."""
    s = np.sqrt(rho * T / M)
    Phi = np.conj(np.swapaxes(PhiH, -1, -2))                # (K, T, M)
    if Hhat.ndim == 2:
        Hhat = Hhat[None]
    XH = s * np.einsum("ktm,lmn->lktn", Phi, Hhat)          # (L, K, T, N)
    R = Y[:, None] - XH
    r2 = np.sum(np.abs(R) ** 2, axis=(2, 3))                # (L, K)
    PR = np.einsum("kmt,lktn->lkmn", PhiH, R)
    pr2 = np.sum(np.abs(PR) ** 2, axis=(2, 3))
    a_e = np.asarray(sig2e) * rho * T / M
    c_e = a_e / (1.0 + a_e)
    c_e = c_e[:, None] if np.ndim(c_e) else c_e
    return -r2 + c_e * pr2


def m_ustm2T(A1, A2, D, c2, lchunk=32):
    """2T-block prior-free struct-ML over K^2 stacked codewords.
    A1, A2: (L, K, M, N) = Phi^H Y_b, Phi^H Y_{b+1}. Returns (L, K, K) with M = 2 rows."""
    L, K = A1.shape[0], A1.shape[1]
    out = np.empty((L, K, K))
    for i in range(0, L, lchunk):
        a1, a2 = A1[i:i + lchunk], A2[i:i + lchunk]
        s11 = np.sum(np.abs(a1[:, :, 0, :]) ** 2, -1)       # (l, K)
        t11 = np.sum(np.abs(a1[:, :, 1, :]) ** 2, -1)
        s22 = np.sum(np.abs(a2[:, :, 0, :]) ** 2, -1)
        t22 = np.sum(np.abs(a2[:, :, 1, :]) ** 2, -1)
        x00 = np.einsum("lkn,ljn->lkj", a1[:, :, 0, :], np.conj(a2[:, :, 0, :]))
        x11 = np.einsum("lkn,ljn->lkj", a1[:, :, 1, :], np.conj(a2[:, :, 1, :]))
        p = 0.5 * (s11[:, :, None] + s22[:, None, :] + 2 * np.real(x00))
        q = 0.5 * (t11[:, :, None] + t22[:, None, :] + 2 * np.real(x11))
        z11 = np.sum(a1[:, :, 0, :] * np.conj(a1[:, :, 1, :]), -1)
        z22 = np.sum(a2[:, :, 0, :] * np.conj(a2[:, :, 1, :]), -1)
        z12 = np.einsum("lkn,ljn->lkj", a1[:, :, 0, :], np.conj(a2[:, :, 1, :]))
        z21 = np.einsum("ljn,lkn->lkj", a2[:, :, 0, :], np.conj(a1[:, :, 1, :]))
        z = 0.5 * (z11[:, :, None] + z22[:, None, :] + z12 + z21)
        mu1, mu2 = C.eig2x2_psd(p, q, z)
        out[i:i + lchunk] = C.log_F(D, c2 * mu1, c2 * mu2)
    return out


# --------------------------------------------------------------------------- trackers
def lam_from_cov(R, N, rank=None):
    """Trace-normalized (sum lam = N) eigen-spectrum of R; if rank is given, keep top-rank only.
    R: (..., N, N) Hermitian PSD. Returns (Q, lam) with lam: (..., d)."""
    ev, V = np.linalg.eigh(R)
    ev = np.maximum(ev, 0.0)
    if rank is not None:
        ev, V = ev[..., -rank:], V[..., -rank:]
    s = np.sum(ev, axis=-1, keepdims=True)
    lam = N * ev / np.maximum(s, 1e-30)
    return V, lam


def debias(R, nblocks, M, rho, T):
    """Remove the LS-estimation noise floor  E[Nt^H Nt] = M^2/(rho T) I_N  per block."""
    N = R.shape[-1]
    return R - nblocks * (M ** 2 / (rho * T)) * np.eye(N)


# --------------------------------------------------------------------------- simulation
def run_traj(H, snr_list, K, L, W0, reffs, dets, seed, T, verbose=True):
    """Returns err[name] = (S, B) error counts (out of L) and meta."""
    B, M, N = H.shape
    r = 2
    rng_c = np.random.default_rng([seed, 777])
    Phi = C.haar_stiefel(rng_c, K, T, M)                     # fixed constellation for everything
    PhiH = C.H_(Phi)

    # --- exact per-block instantaneous subspace and the non-causal long-term subspace
    U_inst = np.stack([np.linalg.svd(H[b], full_matrices=False)[2][:r].conj().T for b in range(B)])
    C_all = np.einsum("bmn,bmk->nk", H.conj(), H)
    U_LT = {rr: topk_eig(C_all, rr) for rr in reffs}
    Ccum = np.cumsum(np.einsum("bmn,bmk->bnk", H.conj(), H), axis=0)     # centered non-causal windows

    def win_cov(s, e):
        return Ccum[e - 1] - (Ccum[s - 1] if s > 0 else 0.0)

    err = {d: np.zeros((len(snr_list), B)) for d in dets}
    for si, sdb in enumerate(snr_list):
        rho = 10 ** (sdb / 10.0)
        a_pf, c_pf = C.snr_consts(rho, T, M, N / r)          # prior-free struct-ML constant
        _, c_2T = C.snr_consts(rho, 2 * T, M, N / r)         # 2T-block
        rng = np.random.default_rng([seed, si, int(zlib.crc32(np.ascontiguousarray(H).tobytes()) % 10 ** 6)])
        j = rng.integers(0, K, (B, L))
        sX = np.sqrt(rho * T / M)

        # decision-directed state (per repetition)
        bufH = np.zeros((L, W0, M, N), complex)              # rolling LS channel estimates
        nfill = 0
        prevA = None                                         # A of block b-1 (for 2T USTM)
        prev_j = None
        for b in range(B):
            W = C.crandn(rng, L, T, N)
            Y = sX * np.einsum("ltm,mn->ltn", Phi[j[b]], H[b]) + W
            A = A_of(PhiH, Y)                                # (L, K, M, N)
            met = {}

            if "glrt" in dets:
                met["glrt"] = m_glrt(A)
            if "struct" in dets:
                met["struct"] = m_struct(A, N, c_pf)
            if "genieU" in dets:
                met["genieU"] = m_cov(A_of(PhiH, Y @ U_inst[b]), _c(np.full(r, N / r), rho, T, M))
            if "genieH" in dets:
                met["genieH"] = m_coh(Y, PhiH, H[b], rho, T, M, 0.0)
            for d in dets:
                mm = re.fullmatch(r"genieLT(\d+)", d)
                if mm:
                    rr = int(mm.group(1))
                    met[d] = m_cov(A_of(PhiH, Y @ U_LT[rr]), _c(np.full(rr, N / rr), rho, T, M))
                mm = re.fullmatch(r"genieLTw(\d+)", d)
                if mm:                       # current long-term subspace, centred NON-CAUSAL window
                    rr = int(mm.group(1))
                    Rw = win_cov(max(0, b - W0 // 2), min(B, max(b + W0 // 2, W0)))
                    Qw = topk_eig(Rw, rr)
                    met[d] = m_cov(A_of(PhiH, Y @ Qw), _c(np.full(rr, N / rr), rho, T, M))
                if d == "genieCovw":         # full second-order knowledge (report only, not D-14)
                    Rw = win_cov(max(0, b - W0 // 2), min(B, max(b + W0 // 2, W0)))
                    Qw, lamw = lam_from_cov(Rw, N)
                    met[d] = m_cov(A_of(PhiH, Y @ Qw), _c(lamw, rho, T, M))

            # ---- oracle-window (causal in H) long-term tracking: upper bound on any causal tracker
            if b >= W0:
                Rw = np.einsum("bmn,bmk->nk", H[b - W0:b].conj(), H[b - W0:b])
                for d in dets:
                    mm = re.fullmatch(r"ltO(\d+)", d)
                    if mm:
                        rr = int(mm.group(1))
                        Q, lam = lam_from_cov(Rw, N, rank=rr)
                        met[d] = m_cov(A_of(PhiH, Y @ Q), _c(lam, rho, T, M))
                    mm = re.fullmatch(r"ltOs(\d+)", d)
                    if mm:
                        rr = int(mm.group(1))
                        Q, _ = lam_from_cov(Rw, N, rank=rr)
                        _, c_e = C.snr_consts(rho, T, M, N / rr)
                        met[d] = m_struct(A_of(PhiH, Y @ Q), rr, c_e)
                if "covO" in dets:
                    Q, lam = lam_from_cov(Rw, N)
                    met["covO"] = m_cov(A_of(PhiH, Y @ Q), _c(lam, rho, T, M))
                if "predO" in dets:
                    num = np.sum([np.vdot(H[t - 1], H[t]) for t in range(b - W0 + 1, b)])
                    den = np.sum([np.linalg.norm(H[t - 1]) ** 2 for t in range(b - W0 + 1, b)])
                    ah = num / max(den, 1e-30)
                    res = np.mean([np.linalg.norm(H[t] - ah * H[t - 1]) ** 2
                                   for t in range(b - W0 + 1, b)]) / (M * N)
                    met["predO"] = m_coh(Y, PhiH, ah * H[b - 1], rho, T, M, float(res))

            # ---- decision-directed tracking (realistic; uses only Y and its own decisions)
            ref_dd = None
            if nfill >= max(W0 // 2, 2):
                Rd = np.einsum("lwmn,lwmk->lnk", bufH[:, :nfill].conj(), bufH[:, :nfill])
                Rd = debias(Rd, nfill, M, rho, T)
                Qd, lamd = lam_from_cov(Rd, N)
                # the tracking loop's own decision (soft conditional ML); computed even if covD
                # is not in --detectors so that every DD variant sees the same channel track
                ref_dd = m_cov(A_of(PhiH, np.einsum("ltn,lnd->ltd", Y, Qd)), _c(lamd, rho, T, M))
                if "covD" in dets:
                    met["covD"] = ref_dd
                for d in dets:
                    mm = re.fullmatch(r"ltD(\d+)", d)
                    if mm:
                        rr = int(mm.group(1))
                        Q, lam = lam_from_cov(Rd, N, rank=rr)
                        Yq = np.einsum("ltn,lnd->ltd", Y, Q)
                        met[d] = m_cov(A_of(PhiH, Yq), _c(lam, rho, T, M))
                if "predD" in dets and nfill >= 2:
                    Hp, Hc = bufH[:, nfill - 2], bufH[:, nfill - 1]
                    num = np.einsum("lmn,lmn->l", Hp.conj(), Hc)
                    den = np.sum(np.abs(Hp) ** 2, axis=(1, 2))
                    ah = num / np.maximum(den, 1e-30)
                    Hhat = ah[:, None, None] * Hc
                    s2 = np.maximum(np.mean(np.abs(Hc - ah[:, None, None] * Hp) ** 2, axis=(1, 2))
                                    - M / (rho * T), 1e-6)
                    met["predD"] = m_coh(Y, PhiH, Hhat, rho, T, M, s2)

            # ---- 2T-block USTM: decide the pair at odd b, charge both blocks
            if "ustm2T" in dets:
                if b % 2 == 1 and prevA is not None:
                    mp = m_ustm2T(prevA, A, N, c_2T)          # (L, K, K)
                    flat = np.argmax(mp.reshape(L, -1), axis=1)
                    k1, k2 = flat // K, flat % K
                    err["ustm2T"][si, b - 1] += np.sum(k1 != prev_j)
                    err["ustm2T"][si, b] += np.sum(k2 != j[b])
                elif b % 2 == 0 and b == B - 1:
                    err["ustm2T"][si, b] += np.sum(np.argmax(m_struct(A, N, c_pf), 1) != j[b])

            # ---- decisions, error counting, and the decision-directed update
            for d, mval in met.items():
                err[d][si, b] += np.sum(np.argmax(mval, axis=1) != j[b])
            ref = ref_dd if ref_dd is not None else m_struct(A, N, c_pf)
            jhat = np.argmax(ref, axis=1)
            Hls = np.sqrt(M / (rho * T)) * A[np.arange(L), jhat]        # LS estimate from own decision
            if nfill < W0:
                bufH[:, nfill] = Hls
                nfill += 1
            else:
                bufH[:, :-1] = bufH[:, 1:]
                bufH[:, -1] = Hls
            prevA, prev_j = A, j[b]
        if verbose:
            top = {d: round(float(err[d][si].sum() / (B * L)), 4) for d in sorted(err)}
            print(f"    snr {sdb:+.0f} dB  SER {top}", flush=True)
    return err


def _c(lam, rho, T, M):
    a = np.asarray(lam) * rho * T / M
    return a / (1.0 + a)


# --------------------------------------------------------------------------- statistics
def ser_curve(errs, L, valid=None):
    """errs: (S, B) counts. valid: block mask (skip warm-up blocks). -> (S,) SER."""
    m = np.ones(errs.shape[1], bool) if valid is None else valid
    return errs[:, m].sum(1) / (L * m.sum())


def boot_cross(errs, L, snr, nboot=400, target=1e-2, seed=0, valid=None):
    """Cluster bootstrap over blocks -> (lo, hi) of the SNR at SER = target."""
    rng = np.random.default_rng(seed)
    idx0 = np.arange(errs.shape[1]) if valid is None else np.where(valid)[0]
    out = []
    for _ in range(nboot):
        idx = rng.choice(idx0, size=len(idx0), replace=True)
        s = errs[:, idx].sum(1) / (L * len(idx))
        x = C.interp_cross(snr, s, target)
        if np.isfinite(x):
            out.append(x)
    if len(out) < nboot // 4:
        return (float("nan"), float("nan"))
    return (float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5)))


# --------------------------------------------------------------------------- self-test
def selftest(T=16, M=2, N=32, r=2, K=64, seed=20260918):
    """T1 (validation against the confirmed gate-(0) numbers): on a synthetic ceiling-model channel
       H_b = G_b U_b^H (U_b Haar per block, G_b iid CN(0, N/r)) the struct-ML / exact-U genie GAP must
       reproduce t6_ceiling's 3.59 dB (16,2,32,2; K=64).  Absolute dB shift by ~0.5 dB with the
       constellation draw, so only the GAP is checked.
       T2: on a CONSTANT channel the 2T-block USTM must beat per-block struct-ML (2T reference)."""
    rng = np.random.default_rng([seed, 4242])
    snr = [-14, -12, -10, -8, -6, -4]
    B, L = 96, 200
    G = np.sqrt(N / r) * C.crandn(rng, B, M, r)
    U = C.haar_stiefel(rng, B, N, r)
    H = normalize(np.einsum("bmr,bnr->bmn", G, U.conj()), "sequence")
    e = run_traj(H, snr, K, L, 16, [4], ["struct", "genieU"], seed, T, verbose=False)
    x = {d: C.interp_cross(snr, ser_curve(e[d], L)) for d in e}
    gap = x["struct"] - x["genieU"]
    ok1 = abs(gap - 3.59) < 0.6
    print(f"[T1] struct {x['struct']:.2f} dB, genieU {x['genieU']:.2f} dB, gap {gap:.2f} dB "
          f"(t6_ceiling: 3.59) -> {'ok' if ok1 else 'FAILED'}")
    snr2 = [-12, -10, -8, -6]
    B, L = 24, 80
    H0 = np.sqrt(N / r) * C.crandn(rng, 1, M, r)
    U0 = C.haar_stiefel(rng, 1, N, r)
    Hc = normalize(np.repeat(np.einsum("bmr,bnr->bmn", H0, U0.conj()), B, 0), "sequence")
    e = run_traj(Hc, snr2, K, L, 8, [4], ["struct", "ustm2T"], seed, T, verbose=False)
    s1, s2 = ser_curve(e["struct"], L), ser_curve(e["ustm2T"], L)
    ok2 = bool(np.all(s2 <= s1 + 1e-9))
    print(f"[T2] constant channel: struct {np.round(s1, 4)} vs ustm2T {np.round(s2, 4)} "
          f"-> {'ok' if ok2 else 'FAILED'}")
    print("selftest ok" if (ok1 and ok2) else "selftest FAILED")
    return ok1 and ok2


# --------------------------------------------------------------------------- driver
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--H", nargs="+", default=["t6_q12_rt_H_q15.json"])
    ap.add_argument("--traj", default="", help="comma-separated substrings; empty = all")
    ap.add_argument("--snr", default="-14:0:2")
    ap.add_argument("--T", type=int, default=16)
    ap.add_argument("--K", type=int, default=64)
    ap.add_argument("--L", type=int, default=200, help="noise repetitions per block")
    ap.add_argument("--B", type=int, default=0, help="blocks to use (0 = all)")
    ap.add_argument("--W", type=int, default=16, help="tracking window (blocks)")
    ap.add_argument("--reffs", default="2,4,8")
    ap.add_argument("--detectors", default="all")
    ap.add_argument("--norm", default="sequence", choices=["sequence", "block"])
    ap.add_argument("--seed", type=int, default=20260918)
    ap.add_argument("--nboot", type=int, default=400)
    ap.add_argument("--out", default="t6_q15.json")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        sys.exit(0 if selftest() else 1)

    if args.quick:
        args.snr, args.L, args.B, args.nboot = "-10:-6:4", 16, 24, 60
    snr = C.parse_snr(args.snr)
    dets = DEFAULT_SET.split(",") if args.detectors == "all" else args.detectors.split(",")
    reffs = sorted({int(x) for x in args.reffs.split(",")} |
                   {int(m.group(1)) for d in dets
                    for m in [re.fullmatch(r"genieLTw?(\d+)", d)] if m})

    runs = load_H(args.H)
    names = [k for k in runs if not args.traj or any(t in k for t in args.traj.split(","))]
    res = dict(meta=dict(date="2026-09-18", seed=args.seed, snr_dB=snr, K=args.K, L=args.L,
                         T=args.T, W=args.W, norm=args.norm, detectors=dets, reffs=reffs,
                         numpy=np.__version__, H_files=args.H), traj={})
    t0 = time.time()
    for name in names:
        H = normalize(runs[name], args.norm)
        if args.B:
            H = H[:args.B]
        print(f"[traj] {name}  B={H.shape[0]} M={H.shape[1]} N={H.shape[2]}  ({time.time()-t0:.0f}s)", flush=True)
        err = run_traj(H, snr, args.K, args.L, args.W, reffs, dets, args.seed, args.T)
        valid = np.zeros(H.shape[0], bool)
        valid[args.W:] = True                                # warm-up blocks excluded for ALL detectors
        entry = dict(B=int(H.shape[0]), n_valid=int(valid.sum()), warmup=int(args.W), det={})
        for d in dets:
            s = ser_curve(err[d], args.L, valid)
            x = C.interp_cross(snr, s)
            lo, hi = boot_cross(err[d], args.L, snr, args.nboot, seed=args.seed, valid=valid)
            entry["det"][d] = {"ser": [float(v) for v in s], "snr_at_1e-2": float(x),
                               "ci": [lo, hi], "nerr": [int(v) for v in err[d][:, valid].sum(1)],
                               "err_block": [[int(v) for v in row] for row in err[d]]}
        res["traj"][name] = entry
        with open(args.out, "w") as f:
            json.dump(res, f, indent=1)
    print("saved", args.out, f"({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
