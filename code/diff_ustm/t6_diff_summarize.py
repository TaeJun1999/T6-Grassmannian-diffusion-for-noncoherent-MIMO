#!/usr/bin/env python3
"""Summarize t6_diff_results.json: SER-1e-2 operating points (dB) per scheme, gaps vs per-block struct-ML
and genie, and SER curves. Usage: python3 t6_diff_summarize.py t6_diff_results.json"""
import json, sys, math
import numpy as np

fn = sys.argv[1] if len(sys.argv) > 1 else 't6_diff_results.json'
d = json.load(open(fn))
lines = []
P = lines.append
P(f"# Q-06 differential USTM pre-estimate — {d['meta']['date']}, seed {d['meta']['seed']}, R={d['meta']['R']} runs x B={d['meta']['B']} blocks per SNR, n_pb={d['meta']['n_pb']} (K>256: half)")
P("alpha = block-to-block AR(1) coefficient of G_b (U fixed); Jakes: alpha = J0(2 pi f_D T T_s).")
P("Schemes: struct = per-block structure-aware ML (prior-free); genie = per-block, U known; diff = HS differential receiver on raw Y;")
P("         diffU = differential on Y U (U known); diffS2 = differential on Y U_hat, U_hat = top-r right sing. vecs of [Y_{b-1}; Y_b] (prior-free).")
P("Operating point = SNR (dB) at SER 1e-2 (log-linear interpolation; nan = no crossing on the grid, i.e. error floor above 1e-2). MC precision ~ +-0.3-0.5 dB.")
P("")
for c in d['configs']:
    T, M, N, r = c['T'], c['M'], c['N'], c['r']
    P(f"## (T,M,N,r) = ({T},{M},{N},{r}), gamma = {c['gamma']:.0f}, SNR grid {c['snr_dB'][0]:.0f}..{c['snr_dB'][-1]:.0f} dB")
    for K, e in c['K'].items():
        x = e['snr_at_ser1em2']
        P(f"### K = {K} (rate {math.log2(int(K))/T:.3f} bit/channel use)")
        P(f"{'scheme':<16}{'op. point [dB]':>16}{'vs struct [dB]':>16}{'vs genie [dB]':>16}")
        for nm in x:
            v = x[nm]
            vs_s = v - x['struct'] if not math.isnan(v) else float('nan')
            vs_g = v - x['genie'] if not math.isnan(v) else float('nan')
            P(f"{nm:<16}{v:>16.2f}{vs_s:>16.2f}{vs_g:>16.2f}")
        P("")
        # SER table (compact)
        rows = e['rows']
        keys = [k for k in rows[0] if k.startswith('SER_')]
        P("SNR[dB] " + " ".join(f"{k[4:]:>13}" for k in keys))
        for row in rows:
            P(f"{row['snr_dB']:7.1f} " + " ".join(f"{row[k]:>13.4f}" for k in keys))
        P("")
txt = "\n".join(lines)
open(fn.replace('.json', '_summary.txt'), 'w').write(txt)
print(txt)

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    ncfg = len(d['configs']); nK = len(d['configs'][0]['K'])
    fig, axes = plt.subplots(ncfg, nK, figsize=(6.2 * nK, 4.6 * ncfg), squeeze=False)
    for i, c in enumerate(d['configs']):
        for j, (K, e) in enumerate(c['K'].items()):
            ax = axes[i][j]; rows = e['rows']; snr = [x['snr_dB'] for x in rows]
            style = {'struct': ('k', '-', 'per-block struct-ML'), 'genie': ('k', '--', 'per-block genie-U')}
            for k in [k for k in rows[0] if k.startswith('SER_')]:
                nm = k[4:]
                ser = np.array([x[k] for x in rows]); ser = np.where(ser > 0, ser, np.nan)
                if nm in style:
                    col, ls, lab = style[nm]
                    ax.semilogy(snr, ser, color=col, ls=ls, lw=2, label=lab)
                else:
                    var, al = nm.split('_a')
                    col = {'diff': 'C3', 'diffU': 'C0', 'diffS2': 'C2'}[var]
                    ls = {'1.0': '-', '0.984': '--', '0.9': '-.', '0.5': ':', '0.0': (0, (1, 1))}[al]
                    ax.semilogy(snr, ser, color=col, ls=ls, lw=1.4, label=f"{var} a={al}")
            ax.axhline(1e-2, color='gray', lw=0.8)
            ax.set_ylim(1e-4, 1); ax.set_xlabel('SNR [dB]'); ax.set_ylabel('SER'); ax.grid(True, which='both', alpha=0.3)
            ax.set_title(f"(T,M,N,r)=({c['T']},{c['M']},{c['N']},{c['r']}), K={K}")
            ax.legend(fontsize=7, ncol=2)
    fig.tight_layout(); fig.savefig(fn.replace('.json', '_plots.png'), dpi=130)
    print("plot saved")
except Exception as ex:
    print("plot skipped:", ex)
