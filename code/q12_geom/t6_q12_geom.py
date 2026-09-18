#!/usr/bin/env python3
"""
t6_q12_geom.py -- Q-12 (a): geometric pre-check of the D-11 regime
  "path gains G_b decorrelate across blocks, receive-side subspace S_b persists for tens of blocks".

Model [our construction; physics = static scatterers, moving UE, one-ring-type local scattering]
  2-D geometry, BS ULA (N elements, d = lambda/2) at the origin along the x-axis, a_rx(theta)_n = exp(j pi n cos theta).
  UE (M-element ULA along x, a_tx(psi)_m = exp(j pi m cos psi)) at p(t) = p0 + v t (cos phi_v, sin phi_v).
  Cluster 1 (direct, single bounce): static scatterers q_s ~ uniform in a strip around the trajectory;
     sub-ray s: AoA_BS = angle(q_s) [fixed], path length L_s(t) = |q_s| + |q_s - p(t)|,
     weight w_s(t) = exp(-|q_s - p(t)|^2 / (2 rho_loc^2))  (local-scattering window that moves with the UE,
     scatterers themselves static -> Doppler spread from the UE-scatterer leg).
  Cluster 2 (far cluster c2 seen through the same local scatterers, double bounce):
     AoA_BS = angle(c2) + small spread delta_s ~ N(0, sigma_theta2^2) [fixed],
     L_s(t) = |c2| + |c2 - q_s| + |q_s - p(t)|, same weights; power P2 relative to cluster 1.
  Sub-ray reflection coefficients beta_s ~ CN(0,1) fixed per run; phase = -2 pi L_s(t) / lambda.
  H_b = sum_c sum_s g_{c,s}(t_b) a_tx(psi_s(t_b)) a_rx(theta_{c,s})^T   (M x N), blocks at t_b = b Delta.

Metrics (per block spacing Delta, lag k, averaged over blocks and runs)
  subspace: rho_S(k) = ||U_b^H U_{b+k}||_F^2 / r,  U_b = top-r right singular basis of H_b (r = 2)
            d_c(k) = sqrt(1 - rho_S(k)) (normalised chordal distance)
  gain:     G_b = H_b U_b,  G~_{b+k} = H_{b+k} U_b (projected on the CURRENT subspace),
            rho_G(k) = |sum tr(G_b^H G~_{b+k})| / sqrt(sum ||G_b||^2 sum ||G~||^2)  -- complex correlation magnitude,
            i.e. the |alpha| of the AR(1) gain model in T6_diff_ustm_note.md; also the per-cluster analytic
            Jakes reference |J0(2 pi v Delta k / lambda)| (isotropic local scattering).
  rank:     energy fraction of the top-2 eigenvectors of the window covariance sum_b H_b^H H_b (Q-03 check).

Usage
  python3 t6_q12_geom.py --selftest
  python3 t6_q12_geom.py --out t6_q12_geom.json      (~ a few minutes)
"""
import argparse, json, sys, time
import numpy as np
from scipy.special import j0

SEED = 20260917
C0 = 299792458.0


def steer(theta, n):
    """ULA steering vector, half-wavelength spacing, angle from array axis. theta: array (...,), returns (..., n)."""
    return np.exp(1j * np.pi * np.arange(n)[None, :] * np.cos(theta)[:, None])


class Geometry:
    def __init__(self, rng, N=32, M=2, fc=28e9, R=50.0, ue_angle_deg=90.0, phi_v_deg=90.0, v=3.0,
                 rho_loc=1.5, n_scat=60, c2=(-40.0, 30.0), P2=0.5, sigma_theta2_deg=1.0, strip_len=40.0, K_los=0.0):
        self.N, self.M, self.lam, self.v = N, M, C0 / fc, v
        self.rho_loc = rho_loc
        self.K_los = K_los  # specular LOS power relative to the mean power of cluster 1
        ang = np.deg2rad(ue_angle_deg)
        self.p0 = R * np.array([np.cos(ang), np.sin(ang)])
        pv = np.deg2rad(phi_v_deg)
        self.dirv = np.array([np.cos(pv), np.sin(pv)])
        # static scatterers in a strip along the trajectory: length strip_len ahead, width 4 rho_loc
        # density: n_scat within a rho_loc-disc -> total = n_scat * (strip area)/(pi rho_loc^2)
        area = strip_len * 4 * rho_loc
        n_tot = int(np.ceil(n_scat * area / (np.pi * rho_loc ** 2)))
        along = rng.uniform(-2 * rho_loc, strip_len, n_tot)
        across = rng.uniform(-2 * rho_loc, 2 * rho_loc, n_tot)
        perp = np.array([-self.dirv[1], self.dirv[0]])
        self.q = self.p0[None, :] + along[:, None] * self.dirv[None, :] + across[:, None] * perp[None, :]
        self.beta1 = (rng.standard_normal(n_tot) + 1j * rng.standard_normal(n_tot)) / np.sqrt(2)
        self.beta2 = (rng.standard_normal(n_tot) + 1j * rng.standard_normal(n_tot)) / np.sqrt(2)
        self.c2 = np.array(c2, dtype=float)
        self.P2 = P2
        self.theta1 = np.arctan2(self.q[:, 1], self.q[:, 0])                      # fixed AoA per scatterer
        self.theta2 = np.arctan2(self.c2[1], self.c2[0]) + np.deg2rad(sigma_theta2_deg) * rng.standard_normal(n_tot)
        self.L1_bs = np.linalg.norm(self.q, axis=1)                                # BS -> q_s
        self.L2_bs = np.linalg.norm(self.c2) + np.linalg.norm(self.q - self.c2[None, :], axis=1)  # BS -> c2 -> q_s
        self.A1 = steer(self.theta1, N)                                           # (n_tot, N)
        self.A2 = steer(self.theta2, N)

    def channel(self, t):
        """H(t) in C^{M x N} (unnormalised)."""
        p = self.p0 + self.v * t * self.dirv
        dq = self.q - p[None, :]                       # UE -> scatterer
        dist = np.linalg.norm(dq, axis=1)
        w = np.exp(-dist ** 2 / (2 * self.rho_loc ** 2))
        psi = np.arctan2(dq[:, 1], dq[:, 0])           # departure angle at the UE (array along x)
        At = steer(psi, self.M)                        # (n_tot, M)
        g1 = w * self.beta1 * np.exp(-2j * np.pi * (self.L1_bs + dist) / self.lam)
        g2 = np.sqrt(self.P2) * w * self.beta2 * np.exp(-2j * np.pi * (self.L2_bs + dist) / self.lam)
        H = (At * g1[:, None]).T @ self.A1 + (At * g2[:, None]).T @ self.A2   # (M, N)
        if self.K_los > 0:
            # specular LOS ray: AoA = direction of the UE, departure = direction of the BS from the UE
            th = np.arctan2(p[1], p[0]); ps = np.arctan2(-p[1], -p[0]); L = np.linalg.norm(p)
            amp = np.sqrt(self.K_los * np.sum(w ** 2))          # cluster-1 mean power = sum_s w_s^2 (E|beta|^2 = 1)
            H = H + amp * np.exp(-2j * np.pi * L / self.lam) * np.outer(steer(np.array([ps]), self.M)[0], steer(np.array([th]), self.N)[0])
        return H


def run_case(rng, geo_kwargs, deltas, B=48, r=2, n_runs=40, lags=(1, 2, 3, 5, 8, 12, 16, 24, 32, 47)):
    """Returns dict Delta -> metrics (rho_S, d_c, rho_G per lag; top-2 energy fraction)."""
    out = {}
    lags = [k for k in lags if k < B]
    for Delta in deltas:
        num_S = {k: 0.0 for k in lags}; cnt_S = {k: 0 for k in lags}
        num_G = {k: 0.0 + 0j for k in lags}; den_G1 = {k: 0.0 for k in lags}; den_G2 = {k: 0.0 for k in lags}
        rhoS_first = {k: 0.0 for k in lags}  # correlation w.r.t. block 0 only (persistence from a fixed reference)
        top2 = 0.0
        for run in range(n_runs):
            geo = Geometry(rng, **geo_kwargs)
            H = np.stack([geo.channel(b * Delta) for b in range(B)])       # (B, M, N)
            # normalise per block (metrics are scale-free anyway)
            H /= (np.linalg.norm(H, axis=(1, 2), keepdims=True) + 1e-30)
            U = np.empty((B, geo.N, r), complex)
            for b in range(B):
                _, _, Vh = np.linalg.svd(H[b], full_matrices=False)
                U[b] = Vh[:r].conj().T                                         # basis of span(H_b^H)
            R = np.einsum('bmn,bmk->nk', H.conj(), H)                          # sum_b H_b^H H_b
            ev = np.linalg.eigvalsh(R)[::-1]
            top2 += ev[:2].sum() / ev.sum()
            for k in lags:
                for b in range(B - k):
                    ov = np.linalg.norm(U[b].conj().T @ U[b + k]) ** 2 / r
                    num_S[k] += ov; cnt_S[k] += 1
                    G0 = H[b] @ U[b]; Gk = H[b + k] @ U[b]
                    num_G[k] += np.vdot(G0, Gk)
                    den_G1[k] += np.linalg.norm(G0) ** 2; den_G2[k] += np.linalg.norm(Gk) ** 2
                rhoS_first[k] += np.linalg.norm(U[0].conj().T @ U[k]) ** 2 / r
        d = geo_kwargs.get('v', 3.0) * Delta / geo.lam
        out[str(Delta)] = {
            'Delta_s': Delta, 'dist_per_block_lambda': d,
            'rho_S': {str(k): num_S[k] / cnt_S[k] for k in lags},
            'rho_S_from_block0': {str(k): rhoS_first[k] / n_runs for k in lags},
            'd_c': {str(k): float(np.sqrt(max(1 - num_S[k] / cnt_S[k], 0))) for k in lags},
            'rho_G': {str(k): float(abs(num_G[k]) / np.sqrt(den_G1[k] * den_G2[k])) for k in lags},
            'jakes_ref': {str(k): float(abs(j0(2 * np.pi * d * k))) for k in lags},
            'top2_energy_fraction_window': top2 / n_runs,
        }
        print(f"  Delta={Delta*1e3:7.2f} ms (d/lambda={d:6.3f}): rho_G(1)={out[str(Delta)]['rho_G']['1']:.3f} "
              f"[J0 {out[str(Delta)]['jakes_ref']['1']:.3f}]  rho_S(1)={out[str(Delta)]['rho_S']['1']:.4f}  "
              f"rho_S(16)={out[str(Delta)]['rho_S'].get('16', float('nan')):.4f}  top2={out[str(Delta)]['top2_energy_fraction_window']:.4f}",
              flush=True)
    return out


def selftest():
    rng = np.random.default_rng(SEED)
    ok = True
    # T1: single isotropic local-scattering cluster, far BS: gain correlation vs Jakes J0 at small displacement
    geo = Geometry(rng, N=8, R=200.0, rho_loc=1.0, n_scat=400, P2=0.0, strip_len=4.0)
    lam = geo.lam
    ds = [0.05, 0.1, 0.2, 0.3]
    err = []
    for dl in ds:
        acc = 0j; d1 = 0.0; d2 = 0.0
        for trial in range(200):
            g2 = Geometry(rng, N=8, R=200.0, rho_loc=1.0, n_scat=400, P2=0.0, strip_len=4.0,
                          phi_v_deg=rng.uniform(0, 360))
            H0 = g2.channel(0.0); H1 = g2.channel(dl * lam / g2.v)
            acc += np.vdot(H0, H1); d1 += np.linalg.norm(H0) ** 2; d2 += np.linalg.norm(H1) ** 2
        rho = abs(acc) / np.sqrt(d1 * d2)
        ref = abs(j0(2 * np.pi * dl))
        err.append(abs(rho - ref))
        print(f"[T1] d/lambda={dl}: |rho_H|={rho:.3f}  |J0|={ref:.3f}")
    ok &= max(err) < 0.08
    # T2: subspace overlap of two steering vectors vs Fejer kernel
    N = 32
    for du in [0.0, 0.01, 0.03, 0.0625]:
        th0 = np.pi / 2; th1 = np.arccos(np.cos(th0) - du)
        a0, a1 = steer(np.array([th0]), N)[0], steer(np.array([th1]), N)[0]
        ov = abs(np.vdot(a0, a1)) ** 2 / N ** 2
        ref = (np.sin(np.pi * N * du / 2) / (N * np.sin(np.pi * du / 2))) ** 2 if du > 0 else 1.0
        print(f"[T2] du={du}: overlap={ov:.4f} fejer={ref:.4f}")
        ok &= abs(ov - ref) < 1e-9
    # T3: H_b has rank <= M and the projected gain reproduces H_b exactly
    geo = Geometry(rng)
    H = geo.channel(0.0); _, s, Vh = np.linalg.svd(H); U = Vh[:2].conj().T
    rec = (H @ U) @ U.conj().T
    e = np.linalg.norm(H - rec) / np.linalg.norm(H)
    print(f"[T3] rank-2 reconstruction error {e:.2e}")
    ok &= e < 1e-10
    print("selftest ok" if ok else "selftest FAILED")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--out', default='t6_q12_geom.json')
    ap.add_argument('--n_runs', type=int, default=40)
    ap.add_argument('--B', type=int, default=48)
    args = ap.parse_args()
    if args.selftest:
        sys.exit(0 if selftest() else 1)
    rng = np.random.default_rng(SEED)
    deltas = [0.5e-3, 1e-3, 2e-3, 5e-3, 10e-3, 20e-3, 50e-3, 100e-3]
    cases = {
        'main_N32_R50_v3_phi90': dict(N=32, R=50.0, v=3.0, phi_v_deg=90.0),      # tangential motion (max AoA drift)
        'main_N32_R50_v3_phi0': dict(N=32, R=50.0, v=3.0, phi_v_deg=0.0),        # radial-ish motion (min AoA drift)
        'N32_R20_v3_phi90': dict(N=32, R=20.0, v=3.0, phi_v_deg=90.0),
        'N32_R100_v3_phi90': dict(N=32, R=100.0, v=3.0, phi_v_deg=90.0),
        'N16_R50_v3_phi90': dict(N=16, R=50.0, v=3.0, phi_v_deg=90.0),
        'N64_R50_v3_phi90': dict(N=64, R=50.0, v=3.0, phi_v_deg=90.0),
        'main_N32_R50_v3_phi90_rho3': dict(N=32, R=50.0, v=3.0, phi_v_deg=90.0, rho_loc=3.0),  # wider local scattering
    }
    res = {'seed': SEED, 'B': args.B, 'n_runs': args.n_runs, 'deltas_s': deltas, 'cases': {}}
    t0 = time.time()
    for name, kw in cases.items():
        print(f"== {name} {kw}", flush=True)
        res['cases'][name] = {'kwargs': kw, 'results': run_case(rng, kw, deltas, B=args.B, n_runs=args.n_runs)}
        with open(args.out, 'w') as fh:
            json.dump(res, fh, indent=1)
        print(f"   ({time.time()-t0:.0f} s)", flush=True)
    print('wrote', args.out)


if __name__ == '__main__':
    main()
