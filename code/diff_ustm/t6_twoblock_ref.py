#!/usr/bin/env python3
"""Reference for Q-06: if the channel is constant over TWO blocks, the best prior-free use of both blocks (random codebook)
is per-block USTM with block length 2T and K^2 points (same rate log2 K / T). Reports struct-ML and genie-U SER at 2T."""
import sys, json, time, argparse
import numpy as np
sys.path.insert(0, '/mnt/project'); sys.path.insert(0, '/home/claude')
import t6_ceiling as C, t6_diff_ustm as D
ap = argparse.ArgumentParser()
ap.add_argument('--cases', default='16,2,16,2,4096,-8:2:1,2000;32,2,32,2,4096,-14:-4:1,1000')
ap.add_argument('--seed', type=int, default=20260916)
ap.add_argument('--out', default='t6_twoblock_ref.json')
args = ap.parse_args()
res = dict(meta=dict(date='2026-09-16', seed=args.seed, cases=args.cases), cases=[])
t0 = time.time()
for ci, case in enumerate(args.cases.split(';')):
    T, M, N, r, K, spec, n = case.split(',')
    T, M, N, r, K, n = int(T), int(M), int(N), int(r), int(K), int(n)
    gamma = N / r
    snr = C.parse_snr(spec)
    rng = np.random.default_rng([args.seed, 77, ci])
    Phic = C.haar_stiefel(rng, K, T, M); PhicH = C.H_(Phic)
    chunk = int(np.clip(2 ** 23 // (K * N * 2), 1, 128))
    rows = []
    for sdb in snr:
        o = D.sim_perblock(C, T, M, N, r, K, 10 ** (sdb / 10), gamma, n, rng, Phic, PhicH, chunk)
        row = dict(snr_dB=float(sdb), **{f'SER_{k}': e / nn for k, (e, nn) in o.items()}, n=n)
        rows.append(row)
        print(f"   T={T} N={N} K={K} SNR={sdb:5.1f} struct={row['SER_struct']:.4f} genie={row['SER_genie']:.4f} ({time.time()-t0:.0f}s)", flush=True)
    cross = {k: C.interp_cross([x['snr_dB'] for x in rows], [x['SER_' + k] for x in rows]) for k in ('struct', 'genie')}
    print(f"   crossings: {cross}", flush=True)
    res['cases'].append(dict(T=T, M=M, N=N, r=r, K=K, gamma=gamma, n=n, rows=rows, snr_at_ser1em2=cross))
    json.dump(res, open(args.out, 'w'), indent=1)
print("saved", args.out, f"({time.time()-t0:.0f}s)")
