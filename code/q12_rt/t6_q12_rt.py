#!/usr/bin/env python3
"""
Q-12 (b): ray-tracing (Sionna RT 2.1.0, Mitsuba 3.9 LLVM) check of the D-11 regime on the built-in 'munich' scene.
BS: documented position [8.5, 21, 27] (m), ULA N=32 (0.5 lambda, iso, V-pol). UE: ULA M=2, height 1.5 m, straight trajectory
at v = 3 m/s along the locally free street direction (probed). 28 GHz. Blocks at t_b = b Delta; per block the narrowband channel
H_b = sum_p a_p (M x N) from the CIR (los + specular + diffuse, max_depth 3, refraction off, synthetic array, fixed seed so that
diffuse-path phases are spatially consistent -- verified: corr 1.000 at the same position, 0.997 at lambda/10).
Scattering coefficient S of all radio materials is a knob: S=0 (specular only), 0.4, 0.8.
Output: t6_q12_rt.json (H sequences + metrics from t6_q12_metrics.metrics).
"""
import json, time, sys, numpy as np
import sionna.rt as rt
from sionna.rt import load_scene, PlanarArray, Transmitter, Receiver, PathSolver
from t6_q12_metrics import metrics

FC = 28e9; LAM = 3e8 / FC; V = 3.0; BS = [8.5, 21.0, 27.0]; SEED = 20260917
scene = load_scene(rt.scene.munich); scene.frequency = FC
scene.tx_array = PlanarArray(num_rows=1, num_cols=2, vertical_spacing=0.5, horizontal_spacing=0.5, pattern="iso", polarization="V")
scene.rx_array = PlanarArray(num_rows=1, num_cols=32, vertical_spacing=0.5, horizontal_spacing=0.5, pattern="iso", polarization="V")
tx = Transmitter(name="ue", position=[0.0, 0.0, 1.5]); rx = Receiver(name="bs", position=BS); scene.add(tx); scene.add(rx)
solver = PathSolver()

def set_S(S):
    for mat in scene.radio_materials.values(): mat.scattering_coefficient = float(S)

def H_at(pos, diffuse, samples=3 * 10 ** 5, max_depth=3):
    tx.position = [float(pos[0]), float(pos[1]), 1.5]
    p = solver(scene, max_depth=max_depth, los=True, specular_reflection=True, diffuse_reflection=diffuse, refraction=False,
               synthetic_array=True, samples_per_src=samples, seed=1)
    a, tau = p.cir(normalize_delays=False, out_type='numpy')
    if a.shape[-2] == 0: return np.zeros((2, 32), complex), 0, False
    H = a[0, :, 0, :, :, 0].sum(-1).T
    d = np.linalg.norm(np.array(BS) - np.array([pos[0], pos[1], 1.5]))
    los = bool(np.any(np.abs(tau.reshape(-1) - d / 3e8) < 1e-9))
    return H, int(a.shape[-2]), los

def probe_direction(start):
    """Pick the axis direction along which the specular-only channel power stays within 6 dB over 3 m (street direction)."""
    set_S(0.0); best = None
    for d in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
        pws = []
        for s in np.linspace(0, 3.0, 6):
            H, n, _ = H_at((start[0] + d[0] * s, start[1] + d[1] * s), False, samples=5 * 10 ** 4, max_depth=2)
            pws.append(np.linalg.norm(H) ** 2)
        pws = np.array(pws); ok = pws.min() > 0
        spread = 10 * np.log10(pws.max() / max(pws.min(), 1e-30)) if ok else 99
        if best is None or spread < best[1]: best = (d, spread)
    return best

def run(name, start, direction, Delta, B, S):
    set_S(S); diffuse = S > 0
    H = []; npaths = []; loss = []
    t0 = time.time()
    for b in range(B):
        pos = (start[0] + direction[0] * V * b * Delta, start[1] + direction[1] * V * b * Delta)
        Hb, n, los = H_at(pos, diffuse); H.append(Hb); npaths.append(n); loss.append(los)
    H = np.stack(H)
    m = metrics(H, W=16, lags=(1, 2, 4, 8, 16, 32, 47))
    # specular-only reference power fraction at the start point
    set_S(0.0); Hs, _, _ = H_at(start, False); set_S(S)
    m.update({'name': name, 'start': list(start), 'direction': list(direction), 'Delta_s': Delta, 'S': S, 'B': B,
              'dist_per_block_lambda': V * Delta / LAM, 'los_fraction_blocks': float(np.mean(loss)), 'npaths_mean': float(np.mean(npaths)),
              'specular_power_fraction': float(np.linalg.norm(Hs) ** 2 / max(np.linalg.norm(H[0]) ** 2, 1e-30)),
              'mean_power': float(np.mean(np.linalg.norm(H, axis=(1, 2)) ** 2)), 'time_s': time.time() - t0,
              'H_re': H.real.tolist(), 'H_im': H.imag.tolist()})
    print(f"{name} S={S} D={Delta*1e3:g}ms B={B}: d/lam={m['dist_per_block_lambda']:.2f} los={m['los_fraction_blocks']:.2f} spec_frac={m['specular_power_fraction']:.2f} "
          f"rhoG(1)={m['rho_G']['1']:.3f} rhoG(4)={m['rho_G']['4']:.3f} rhoS(1)={m['rho_S']['1']:.3f} r95={m['r_eff95']} PEF2(1)={m['PEF_r2'].get('1',float('nan')):.3f} "
          f"PEF2(32)={m['PEF_r2'].get('32',float('nan')):.3f} PEFe(32)={m['PEF_reff95'].get('32',float('nan')):.3f} LT(32)={m['LT_overlap_reff95'].get('32',float('nan')):.3f} ({m['time_s']:.0f}s)", flush=True)
    return m

if __name__ == '__main__':
    points = {'LOS_A': (40.0, 48.0), 'NLOS_A': (-8.0, 48.0), 'NLOS_B': (8.0, 40.0)}
    out = {'scene': 'munich', 'BS': BS, 'fc': FC, 'v': V, 'seed': 1, 'runs': []}
    dirs = {}
    for name, st in points.items():
        d, spread = probe_direction(st); dirs[name] = d; print(name, 'direction', d, 'power spread dB', round(spread, 1), flush=True)
    plan = []
    for name in points:
        for S in (0.0, 0.4, 0.8): plan.append((name, 10e-3, 80, S))
        plan.append((name, 2e-3, 80, 0.4)); plan.append((name, 100e-3, 48, 0.4))
    for (name, Delta, B, S) in plan:
        out['runs'].append(run(name, points[name], dirs[name], Delta, B, S))
        json.dump(out, open('t6_q12_rt.json', 'w'))
    print('done')
