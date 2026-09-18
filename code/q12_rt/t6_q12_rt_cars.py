#!/usr/bin/env python3
"""Q-12 (b) part 2: Munich + car clutter (24 low-poly cars, metal, S=0.4) placed within +-3 m of the UE trajectory.
Same BS/arrays/solver settings as t6_q12_rt.py. Output t6_q12_rt_cars.json."""
import json, time, os, numpy as np
import sionna.rt as rt, mitsuba as mi
from sionna.rt import load_scene, PlanarArray, Transmitter, Receiver, PathSolver, SceneObject, RadioMaterial
from t6_q12_metrics import metrics
FC=28e9; LAM=3e8/FC; V=3.0; BS=[8.5,21.0,27.0]
scene = load_scene(rt.scene.munich); scene.frequency = FC
scene.tx_array = PlanarArray(num_rows=1, num_cols=2, vertical_spacing=0.5, horizontal_spacing=0.5, pattern="iso", polarization="V")
scene.rx_array = PlanarArray(num_rows=1, num_cols=32, vertical_spacing=0.5, horizontal_spacing=0.5, pattern="iso", polarization="V")
for m in scene.radio_materials.values(): m.scattering_coefficient = 0.4
metal = RadioMaterial(name="car_metal", thickness=0.01, relative_permittivity=1.0, conductivity=1e7, scattering_coefficient=0.4); scene.add(metal)
car_ply = os.path.join(os.path.dirname(rt.__file__), 'scenes', 'low_poly_car.ply')
tx = Transmitter(name="ue", position=[0.0,0.0,1.5]); rx = Receiver(name="bs", position=BS); scene.add(tx); scene.add(rx)
solver = PathSolver(); cars = []
def place_cars(start, d, n=24, seed=1):
    global cars
    if cars: scene.edit(remove=[c for c in cars]); cars = []
    rng = np.random.default_rng(seed); perp = np.array([-d[1], d[0]]); start = np.array(start); d = np.array(d, float)
    objs = []
    for i in range(n):
        s = rng.uniform(-4, 10); lat = rng.choice([-2.5, 2.5]) + rng.uniform(-0.5, 0.5); p = start + d * s + perp * lat
        o = SceneObject(fname=car_ply, name=f"car{i}", radio_material=metal); objs.append((o, p, rng.uniform(0, 2*np.pi)))
    scene.edit(add=[o for o, _, _ in objs])
    for o, p, yaw in objs:
        o.position = mi.Point3f(float(p[0]), float(p[1]), 0.8); o.orientation = mi.Point3f(float(yaw), 0.0, 0.0)
    cars = [o for o, _, _ in objs]
def H_at(pos, diffuse=True):
    tx.position = [float(pos[0]), float(pos[1]), 1.5]
    p = solver(scene, max_depth=3, los=True, specular_reflection=True, diffuse_reflection=diffuse, refraction=False, synthetic_array=True, samples_per_src=3*10**5, seed=1)
    a, tau = p.cir(normalize_delays=False, out_type='numpy')
    if a.shape[-2] == 0: return np.zeros((2, 32), complex), 0, False
    dd = np.linalg.norm(np.array(BS) - np.array([pos[0], pos[1], 1.5]))
    return a[0, :, 0, :, :, 0].sum(-1).T, int(a.shape[-2]), bool(np.any(np.abs(tau.reshape(-1) - dd/3e8) < 1e-9))
def run(name, start, d, Delta, B):
    t0 = time.time(); H = []; los = []
    for b in range(B):
        Hb, n, l = H_at((start[0] + d[0]*V*b*Delta, start[1] + d[1]*V*b*Delta)); H.append(Hb); los.append(l)
    H = np.stack(H); m = metrics(H, W=16, lags=(1, 2, 4, 8, 16, 32, 47))
    Hs, _, _ = H_at(start, False)
    m.update({'name': name, 'start': list(start), 'direction': list(d), 'Delta_s': Delta, 'B': B, 'clutter': '24 cars, metal, S=0.4',
              'dist_per_block_lambda': V*Delta/LAM, 'los_fraction_blocks': float(np.mean(los)),
              'specular_power_fraction': float(np.linalg.norm(Hs)**2 / max(np.linalg.norm(H[0])**2, 1e-30)), 'time_s': time.time()-t0,
              'H_re': H.real.tolist(), 'H_im': H.imag.tolist()})
    print(f"{name}+cars D={Delta*1e3:g}ms: d/lam={m['dist_per_block_lambda']:.2f} los={m['los_fraction_blocks']:.2f} spec_frac={m['specular_power_fraction']:.2f} "
          f"rhoG(1)={m['rho_G']['1']:.3f} rhoG(4)={m['rho_G']['4']:.3f} rhoS(1)={m['rho_S']['1']:.3f} r95={m['r_eff95']} PEF2(1)={m['PEF_r2'].get('1',float('nan')):.3f} "
          f"PEF2(32)={m['PEF_r2'].get('32',float('nan')):.3f} PEFe(32)={m['PEF_reff95'].get('32',float('nan')):.3f} LT(32)={m['LT_overlap_reff95'].get('32',float('nan')):.3f} ({m['time_s']:.0f}s)", flush=True)
    return m
if __name__ == '__main__':
    out = {'scene': 'munich + cars', 'BS': BS, 'runs': []}
    for name, st, d in [('NLOS_A', (-8.0, 48.0), (0, -1)), ('LOS_far', (72.0, 72.0), None)]:
        if d is None:
            # probe direction for LOS_far (specular only, quick)
            best = None
            for dd in [(1,0),(-1,0),(0,1),(0,-1)]:
                pws = [np.linalg.norm(H_at((st[0]+dd[0]*s, st[1]+dd[1]*s), False)[0])**2 for s in np.linspace(0,3,4)]
                sp = 10*np.log10(max(pws)/max(min(pws),1e-30)) if min(pws) > 0 else 99
                if best is None or sp < best[1]: best = (dd, sp)
            d = best[0]; print('LOS_far direction', d, round(best[1],1), flush=True)
        place_cars(st, d)
        for Delta, B in ([(10e-3, 80), (2e-3, 80)] if name == 'NLOS_A' else [(10e-3, 80)]):
            out['runs'].append(run(name, st, d, Delta, B)); json.dump(out, open('t6_q12_rt_cars.json', 'w'))
    print('done')
