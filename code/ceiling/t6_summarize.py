#!/usr/bin/env python3
"""Summarize t6_ceiling_results.json: tables in bit/s/Hz, gap fractions, SNR gaps at SER 1e-2, CSV and plots."""
import json
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from t6_ceiling import interp_cross

path = sys.argv[1] if len(sys.argv) > 1 else "t6_ceiling_results.json"
d = json.load(open(path))
lines = []
csv = ["config,input,snr_dB,quantity,value_bps_hz,stderr_bps_hz"]


def P(s=""):
    lines.append(s); print(s)


for cfg in d["configs"]:
    T, M, N, r, gam = cfg["T"], cfg["M"], cfg["N"], cfg["r"], cfg["gamma"]
    f = T * np.log(2.0)
    tag = f"({T},{M},{N},{r})"
    P(f"\n=== config (T,M,N,r) = {tag}, gamma = {gam} ===")
    P("--- continuous isotropic USTM input [bit/s/Hz] ---")
    P(f"{'SNR':>5} {'R_genie':>8} {'GMI_svd':>8} {'gapUB':>7} {'gapUB%':>7} {'I(U;Y|X)':>9} {'C_coh':>7}")
    for row in cfg["ustm"]:
        Rg, Gs, Cc, Iu = row["R_genie"] / f, row["GMI_svd"] / f, row["C_coh"] / f, row["I_U_given_X"] / f
        P(f"{row['snr_dB']:5.0f} {Rg:8.3f} {Gs:8.3f} {Rg-Gs:7.3f} {100*(Rg-Gs)/Rg:6.1f}% {Iu:9.3f} {Cc:7.3f}")
        for qn, qv, se in [("R_genie", row["R_genie"], row["se_genie"]), ("GMI_svd", row["GMI_svd"], row["se_svd"]),
                           ("C_coh", row["C_coh"], row["se_coh"]), ("I_U_given_X", row["I_U_given_X"], row["se_IUX"])]:
            csv.append(f"{tag},USTM,{row['snr_dB']},{qn},{qv/f:.5f},{se/f:.5f}")
    for K, rows in cfg["finite"].items():
        if "snr" in K:
            continue
        P(f"--- finite random-Haar constellation K = {K} (cap {np.log2(int(K))/T:.3f} bit/s/Hz), n = {rows[0]['n']} blocks/point ---")
        P(f"{'SNR':>5} {'genie':>7} {'struct':>7} {'svd':>7} {'glrt':>7} {'iid':>7} | {'gap_str':>7} {'gap%':>6} | SER genie/struct/svd/glrt/iid")
        for row in rows:
            g, s_, v, l, i = (row["I_genie"] / f, row["I_struct"] / f, row["GMI_svd"] / f, row["GMI_glrt"] / f, row["I_iid"] / f)
            gap = g - s_
            P(f"{row['snr_dB']:5.0f} {g:7.3f} {s_:7.3f} {v:7.3f} {l:7.3f} {i:7.3f} | {gap:7.3f} {100*gap/max(g,1e-9):5.1f}% | "
              f"{row['SER_genie']:.4f}/{row['SER_struct']:.4f}/{row['SER_svd']:.4f}/{row['SER_glrt']:.4f}/{row['SER_iid']:.4f}")
            for qn in ["I_genie", "I_struct", "GMI_svd", "GMI_glrt", "I_iid"]:
                sek = "se_" + qn.split("_")[1]
                csv.append(f"{tag},K={K},{row['snr_dB']},{qn},{row[qn]/f:.5f},{row[sek]/f:.5f}")
            for qn in ["SER_genie", "SER_struct", "SER_svd", "SER_glrt", "SER_iid"]:
                csv.append(f"{tag},K={K},{row['snr_dB']},{qn},{row[qn]:.5f},")
        cross = {nm: interp_cross([x["snr_dB"] for x in rows], [x["SER_" + nm] for x in rows])
                 for nm in ["genie", "struct", "svd", "glrt", "iid"]}
        P("SNR @ SER 1e-2 [dB]: " + ", ".join(f"{k}={v:.2f}" for k, v in cross.items()) +
          f"  ->  struct-genie = {cross['struct']-cross['genie']:.2f} dB, svd-genie = {cross['svd']-cross['genie']:.2f} dB, glrt-genie = {cross['glrt']-cross['genie']:.2f} dB")

open("t6_ceiling_summary.txt", "w").write("\n".join(lines) + "\n")
open("t6_ceiling_summary.csv", "w").write("\n".join(csv) + "\n")

# ---- plots
ncfg = len(d["configs"])
fig, axes = plt.subplots(3, ncfg, figsize=(6.2 * ncfg, 12))
axes = np.array(axes).reshape(3, ncfg)
for ci, cfg in enumerate(d["configs"]):
    T = cfg["T"]; f = T * np.log(2.0); tag = f"(T,M,N,r)=({cfg['T']},{cfg['M']},{cfg['N']},{cfg['r']})"
    snr = [x["snr_dB"] for x in cfg["ustm"]]
    ax = axes[0, ci]
    ax.plot(snr, [x["R_genie"] / f for x in cfg["ustm"]], "k-o", ms=3, label="genie: subspace known (exact MI, USTM)")
    ax.plot(snr, [x["GMI_svd"] / f for x in cfg["ustm"]], "b-s", ms=3, label="prior-free SVD-projection (GMI, lower bound)")
    ax.plot(snr, [x["C_coh"] / f for x in cfg["ustm"]], "g--", label="coherent, Gaussian input")
    ax.set_xlabel("SNR $\\rho$ [dB]"); ax.set_ylabel("bit/s/Hz"); ax.set_title(tag + "  continuous USTM"); ax.grid(alpha=.3); ax.legend(fontsize=8)
    ax = axes[1, ci]
    ax.plot(snr, [100 * (x["R_genie"] - x["GMI_svd"]) / x["R_genie"] for x in cfg["ustm"]], "b-s", ms=3,
            label="USTM: $(R_{genie}-GMI_{svd})/R_{genie}$ (upper bound on gap)")
    for K, rows in cfg["finite"].items():
        if "snr" in K:
            continue
        ax.plot([x["snr_dB"] for x in rows], [100 * (x["I_genie"] - x["I_struct"]) / max(x["I_genie"], 1e-9) for x in rows],
                "-o", ms=3, label=f"K={K}: exact gap genie $-$ struct-ML")
    ax.axhline(10, color="r", ls=":", label="gate target 10%")
    ax.set_xlabel("SNR $\\rho$ [dB]"); ax.set_ylabel("rate gap [% of genie]"); ax.set_ylim(0, 100); ax.grid(alpha=.3); ax.legend(fontsize=7)
    ax = axes[2, ci]
    for K, rows in cfg["finite"].items():
        if "snr" in K:
            continue
        for nm, st in [("genie", "k-"), ("struct", "r-"), ("svd", "b--"), ("glrt", "c:"), ("iid", "g-.")]:
            ax.semilogy([x["snr_dB"] for x in rows], [max(x["SER_" + nm], 1e-5) for x in rows], st, label=f"K={K} {nm}")
    ax.axhline(1e-2, color="gray", ls=":")
    ax.set_xlabel("SNR $\\rho$ [dB]"); ax.set_ylabel("uncoded SER"); ax.set_ylim(1e-4, 1); ax.grid(alpha=.3, which="both"); ax.legend(fontsize=6, ncol=2)
plt.tight_layout()
plt.savefig("t6_ceiling_plots.png", dpi=130)
print("\nwrote t6_ceiling_summary.txt, t6_ceiling_summary.csv, t6_ceiling_plots.png")
