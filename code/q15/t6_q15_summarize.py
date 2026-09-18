#!/usr/bin/env python3
"""
T6 Q-15 summarizer (2026-09-18): applies the PRE-REGISTERED reading rule of T6_q15_spec_2026-09-18.md
to the output of t6_q15_detectors.py. The rule is compiled in here on purpose -- a gate that is read
by eye moves after the fact.

Reads one or more t6_q15*.json (e.g. the main run plus the separate ustm2T run), merges them per
trajectory, and reports per trajectory:
    SNR at SER 1e-2 for every detector (recomputed from the per-block error counts, warm-up excluded)
    Headroom = SNR(best REALISTIC baseline)  - SNR(best genieLT r)
    Margin   = SNR(predD)                    - SNR(best genieLT r)
with a PAIRED cluster bootstrap over blocks (the two detectors are resampled on the SAME blocks, so
the CI of the difference is far tighter than combining two marginal CIs).

Groups (fixed):
    realistic baselines : struct, glrt, ltD*, covD, predD, ustm2T      <- what gate (i) must beat
    T6 ceiling          : genieLT*                                     <- best case for ANY LT prior
    upper bounds (report only, never judged) : ltO*, covO, predO, genieU, genieH

Verdict (R1-R4, judged on the NLOS trajectories only -- the low-SNR main regime):
    R1  all NLOS have Headroom CI-upper < 1.0 dB   -> (i) unreachable even with a perfect prior
    R2  predD is the best baseline and Margin <= 0 -> same conclusion, stronger
    R3  some NLOS has Headroom point estimate >= 1.0 dB -> proceed to Phase 2
    R4  0.5 <= Headroom < 1.0                      -> borderline: Q-12 (c) and Q-14 first

Usage:
    python3 t6_q15_summarize.py t6_q15_main.json t6_q15_ustm2T.json --out t6_q15_summary.txt
"""
import argparse, json, os, re, sys
import numpy as np

sys.path[:0] = [os.path.dirname(os.path.abspath(__file__)), os.getcwd(),
                os.path.dirname(os.path.dirname(os.path.abspath(__file__)))]
import t6_ceiling as C      # must be reachable: keep it next to this file or set PYTHONPATH

REALISTIC = ("struct", "glrt", "covD", "predD", "ustm2T")          # plus ltD<r>
CEILING_RE = re.compile(r"genieLTw?(\d+)")
UPPER = ("covO", "predO", "genieU", "genieH", "genieCovw")                       # plus ltO<r>, ltOs<r>
NLOS_KEYS = ("nlos_",)                                              # P5 = the only true NLOS point


def is_realistic(d):
    return d in REALISTIC or re.fullmatch(r"ltD\d+", d) is not None


def is_upper(d):
    return d in UPPER or re.fullmatch(r"ltOs?\d+", d) is not None


def cross_flag(snr, ser, target=1e-2):
    """(value, flag). flag: ok | below (already under target at the lowest SNR) | above (floor)."""
    x = C.interp_cross(snr, ser, target)
    if np.isfinite(x):
        return float(x), "ok"
    if ser[0] < target:
        return float(snr[0]), "below"
    return float(snr[-1]), "above"


def ser_from(err, L, valid):
    return err[:, valid].sum(1) / (L * valid.sum())


def paired_boot(errA, errB, L, snr, valid, nboot=2000, seed=0):
    """CI of SNR(A) - SNR(B) with blocks resampled jointly. Replicates where either detector has no
    crossing are dropped; the drop fraction is returned so a degenerate CI cannot be read as tight."""
    rng = np.random.default_rng(seed)
    idx0 = np.where(valid)[0]
    d = []
    for _ in range(nboot):
        idx = rng.choice(idx0, size=len(idx0), replace=True)
        sA = errA[:, idx].sum(1) / (L * len(idx))
        sB = errB[:, idx].sum(1) / (L * len(idx))
        xA, xB = C.interp_cross(snr, sA), C.interp_cross(snr, sB)
        if np.isfinite(xA) and np.isfinite(xB):
            d.append(xA - xB)
    if len(d) < nboot // 4:
        return float("nan"), float("nan"), 1.0 - len(d) / nboot
    return (float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5)), 1.0 - len(d) / nboot)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("json", nargs="+")
    ap.add_argument("--nboot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260918)
    ap.add_argument("--out", default="t6_q15_summary.txt")
    args = ap.parse_args()

    merged, meta = {}, None
    for p in args.json:
        d = json.load(open(p))
        meta = meta or d["meta"]
        for name, e in d["traj"].items():
            m = merged.setdefault(name, dict(B=e["B"], warmup=e.get("warmup", 16), det={}))
            for dd, v in e["det"].items():
                if "err_block" not in v:
                    print(f"  ! {name}/{dd}: no err_block (old run) -- rerun to get paired CIs")
                    continue
                m["det"][dd] = np.asarray(v["err_block"], float)
    snr = meta["snr_dB"]
    L = meta["L"]
    out = [f"# T6 Q-15 summary (rule: T6_q15_spec_2026-09-18.md §4)",
           f"# runs: {args.json}", f"# seed {meta['seed']}  K {meta['K']}  L {L}  T {meta['T']}  "
           f"W {meta['W']}  norm {meta['norm']}  SNR {snr[0]}..{snr[-1]} dB  numpy {meta['numpy']}", ""]

    verdict = {}
    for name in sorted(merged):
        e = merged[name]
        valid = np.zeros(e["B"], bool)
        valid[e["warmup"]:] = True
        xs, flags = {}, {}
        for d, err in e["det"].items():
            xs[d], flags[d] = cross_flag(snr, ser_from(err, L, valid))
        out.append(f"## {name}   (B={e['B']}, valid={int(valid.sum())} blocks, "
                   f"{int(valid.sum())*L} trials/SNR)")
        out.append(f"{'detector':12s} {'SNR@1e-2':>9s}  {'flag':6s} {'role':10s} "
                   f"{'errors @ crossing-ish':>0s}")
        for d in sorted(xs, key=lambda k: xs[k]):
            role = "baseline" if is_realistic(d) else ("ceiling" if CEILING_RE.fullmatch(d)
                                                       else ("upper" if is_upper(d) else "-"))
            out.append(f"{d:12s} {xs[d]:9.2f}  {flags[d]:6s} {role:10s} "
                       f"{[int(v) for v in e['det'][d][:, valid].sum(1)]}")

        base = {d: xs[d] for d in xs if is_realistic(d) and flags[d] != "above"}
        ceil = {d: xs[d] for d in xs if CEILING_RE.fullmatch(d)}
        if not base or not ceil:
            out.append("  (no judgeable pair)\n")
            continue
        db = min(base, key=base.get)
        dc = min(ceil, key=ceil.get)
        hd = xs[db] - xs[dc]
        lo, hi, drop = paired_boot(e["det"][db], e["det"][dc], L, snr, valid, args.nboot, args.seed)
        out.append(f"  best realistic baseline : {db} ({xs[db]:.2f} dB)")
        out.append(f"  best LT ceiling         : {dc} ({xs[dc]:.2f} dB)")
        out.append(f"  HEADROOM = {hd:+.2f} dB   95% CI [{lo:+.2f}, {hi:+.2f}]  (dropped {drop:.0%})")
        mg = None
        if "predD" in xs and flags["predD"] != "above":
            mg = xs["predD"] - xs[dc]
            mlo, mhi, mdrop = paired_boot(e["det"]["predD"], e["det"][dc], L, snr, valid,
                                          args.nboot, args.seed + 1)
            out.append(f"  MARGIN(predD) = {mg:+.2f} dB   95% CI [{mlo:+.2f}, {mhi:+.2f}]")
        warn = [d for d in (db, dc) if flags[d] != "ok"]
        if warn:
            out.append(f"  ! crossing outside the SNR grid for {warn} -- extend --snr before judging")
        out.append("")
        verdict[name] = dict(headroom=hd, hi=hi, margin=mg, best=db, ceil=dc,
                             nlos=any(k in name for k in NLOS_KEYS), ok=not warn)

    nl = {k: v for k, v in verdict.items() if v["nlos"]}
    out.append("## verdict (NLOS trajectories only; LOS reported but not judged, D-12)")
    if not nl:
        out.append("  no NLOS trajectory in this run -- cannot judge.")
    elif not all(v["ok"] for v in nl.values()):
        out.append("  SNR grid too narrow on at least one judging pair -- rerun before judging.")
    else:
        if all(v["hi"] < 1.0 for v in nl.values()):
            out.append("  R1 FIRES: headroom CI-upper < 1.0 dB on every NLOS trajectory ->")
            out.append("     a perfect long-term-subspace prior cannot deliver gate (i) here.")
            out.append("     Record (i) as MISSED in the main regime; do not start training.")
        elif any(v["headroom"] >= 1.0 for v in nl.values()):
            out.append("  R3 FIRES: headroom >= 1.0 dB on " +
                       ", ".join(k for k, v in nl.items() if v["headroom"] >= 1.0) +
                       " -> proceed to Phase 2; target gain = that headroom.")
        else:
            out.append("  R4: borderline (0.5-1.0 dB) -> Q-12 (c) and Q-14 first, then re-read.")
        if any(v["best"] == "predD" and v["margin"] is not None and v["margin"] <= 0
               for v in nl.values()):
            out.append("  R2 ALSO FIRES: the full-channel predictor beats the LT ceiling ->")
            out.append("     gain prediction substitutes for subspace knowledge in this regime.")
    txt = "\n".join(out)
    print(txt)
    open(args.out, "w").write(txt + "\n")
    print("\nsaved", args.out)


if __name__ == "__main__":
    main()
