#!/usr/bin/env python3
"""
reconstruct_energy.py — the headline test: how well do we reconstruct the neutrino
energy if we just take each channel's average missing energy and add it back to the
deposited energy?

    E_reco = E_dep + (missing-energy correction from the channel averages)

Two granularities:
  --level event   one correction per event, keyed by the event topology (lepton x nP x nPi
                  x nShower). Uses EventTree only -> RUNS ON THE CURRENT FILE, no rerun.
  --level vertex  sum a per-vertex correction over the event's vertices, each keyed by its
                  observable VertexMissTree key. Needs the Run/SubRun branches just added
                  to VertexMissTree (rerun), so EventTree E_dep can be joined per event.

Honesty: a lookup trained and tested on the same events flatters itself. We report the
in-sample number (what was literally asked) AND a 5-fold cross-validated number (each
event corrected by averages built from the OTHER events) — the CV number is the real one.

  python3 reconstruct_energy.py bigsampleV2.root --level event
  python3 reconstruct_energy.py bigsampleV2.root --level vertex --csv reco_channels.csv --plot reco.png
  python3 reconstruct_energy.py bigsampleV2.root --level event --keep-uncontained
"""
import argparse, sys
import numpy as np

try:
    import uproot
except ImportError:
    sys.exit("need uproot:  pip install uproot awkward --break-system-packages")


def find_tree(f, name):
    for k in f.keys():
        if k.split(";")[0].split("/")[-1] == name:
            return f[k]
    return None


def stat_line(frac):
    return (f"bias {frac.mean()*100:+6.2f}%   RMS {frac.std()*100:5.2f}%   "
            f"68% {np.percentile(np.abs(frac),68)*100:5.2f}%")


def kfold_indices(n, k, seed=0):
    rng = np.random.default_rng(seed)
    order = rng.permutation(n)
    return [order[i::k] for i in range(k)]      # k disjoint index sets


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rootfile")
    ap.add_argument("--level", choices=["event", "vertex"], default="event")
    ap.add_argument("--mode", choices=["additive", "multiplicative"], default="additive",
                    help="additive: E_reco=E_dep+<miss>_chan. multiplicative: E_reco=E_dep*k_chan "
                         "(k=1/<E_dep/E_nu>_chan), which scales the correction with energy. "
                         "multiplicative is event-level only.")
    ap.add_argument("--keep-uncontained", action="store_true",
                    help="skip the containment cut (default: contained events only)")
    ap.add_argument("--kfolds", type=int, default=5)
    ap.add_argument("--min-count", type=int, default=1,
                    help="channels with fewer than this many MC entries fall back to the global mean")
    ap.add_argument("--csv", default=None, help="write the channel correction table (mean,std,N)")
    ap.add_argument("--event-csv", default=None,
                    help="dump per-event E_true/E_dep/E_reco + fractional errors (CV), for plot_resolution.C")
    ap.add_argument("--plot", default="reconstruct_energy.png")
    ap.add_argument("--decay-correct", action="store_true",
                    help="use the decay-corrected deposit E_dep' = E_dep + E_decay_EndE - E_decay_dDep "
                         "(recovers invisible decay-neutrino energy). Needs the decay branches (rerun).")
    ap.add_argument("--decay-rest-mass", action="store_true",
                    help="with --decay-correct, use the rest-mass variant (E_decay_mass) instead of EndE")
    ap.add_argument("--decay-nu-inside", action="store_true",
                    help="RECOMMENDED decay recovery: add E_neutrino_inside (true decay-neutrino energy). "
                         "Disjoint from the escape term (escape excludes neutrinos), so no double-count. "
                         "Needs only the existing E_neutrino_inside branch (no rerun).")
    ap.add_argument("--max-escapefrac", type=float, default=None, metavar="F",
                    help="Drop events whose truth leakage fraction EscapeFrac = "
                         "(E_escape_charged + E_escape_neutral)/E_nu exceeds F (e.g. 0.3). "
                         "Removes genuine high-leakage events (ultra-high-E flux-tail events that "
                         "leak most of their energy to escaping neutrons) that the EdgeE containment "
                         "cut misses because neutral leakage deposits nothing near the edges. These "
                         "are not reconstructible and poison the channel means. Needs the existing "
                         "EscapeFrac branch (no rerun).")
    ap.add_argument("--keep-neutron-incoming", action="store_true",
                    help="Keep vertices whose incoming class is a neutron (In_class=='n'). "
                         "By default these are DROPPED: a neutron leaves no track, so an "
                         "n->... vertex (especially invisible n->n scatters) can't be identified "
                         "in real data — keying a correction on it uses truth the detector lacks. "
                         "Vertex level only.")
    ap.add_argument("--neutron-conservation", action="store_true",
                    help="Vertex additive only. At observable p->/pi-> vertices, replace the channel-mean "
                         "correction with the per-vertex conservation estimate "
                         "In_KE - E_binding - OutChargedKE (minus OutShowerKE if present), clamped >=0. "
                         "nu-> and other vertices keep the channel mean. Needs In_KE/E_binding/"
                         "OutChargedKE on VertexMissTree.")
    args = ap.parse_args()
    if args.neutron_conservation and (args.level != "vertex" or args.mode != "additive"):
        sys.exit("--neutron-conservation requires --level vertex --mode additive")
    if args.mode == "multiplicative" and args.level == "vertex":
        sys.exit("multiplicative mode is event-level only (you can't scale a per-vertex E_dep). "
                 "Use --level event --mode multiplicative, or --level vertex --mode additive.")

    f = uproot.open(args.rootfile)

    # ---- load the per-event energy (both levels need E_dep and E_nu per event) ----
    et = find_tree(f, "EventTree")
    if et is None:
        sys.exit("no EventTree")
    need_ev = ["Run", "SubRun", "Event", "E_dep_total", "Gen_nu_E", "CCNC"]
    for b in need_ev:
        if b not in et.keys():
            sys.exit(f"EventTree missing '{b}'")
    has_cont = "Contained" in et.keys()
    has_escf = "EscapeFrac" in et.keys()
    ev = et.arrays(need_ev + (["Contained"] if has_cont else [])
                   + (["EscapeFrac"] if has_escf else [])
                   + [b for b in ["Topo_lepton", "Topo_nP", "Topo_nPi", "Topo_nShower"] if b in et.keys()],
                   library="np")
    ev = {k: np.asarray(v) for k, v in ev.items()}

    enu_all  = ev["Gen_nu_E"]
    edep_raw = np.array(ev["E_dep_total"], dtype=float)   # raw calorimetric deposit = the baseline
    decay_add = np.zeros_like(edep_raw)                   # decay recovery is PART OF THE CORRECTION, not the baseline
    if args.decay_nu_inside:
        if "E_neutrino_inside" not in et.keys():
            sys.exit("--decay-nu-inside needs the E_neutrino_inside branch on EventTree.")
        decay_add = np.asarray(et.arrays(["E_neutrino_inside"], library="np")["E_neutrino_inside"], dtype=float)
        print(f"[decay recovery = E_neutrino_inside: mean {decay_add.mean()*1000:+.1f} MeV/event "
              f"(disjoint from escape — no double-count)]")
    elif args.decay_correct:
        add_br = "E_decay_mass" if args.decay_rest_mass else "E_decay_EndE"
        need_d = [add_br, "E_decay_dDep"]
        miss_d = [b for b in need_d if b not in et.keys()]
        if miss_d:
            sys.exit(f"--decay-correct needs {miss_d} on EventTree — rebuild the module and rerun.")
        has_esc = "E_decay_dEsc" in et.keys()      # subtree escape -> makes EndE recovery neutrino-only
        dd = et.arrays(need_d + (["E_decay_dEsc"] if has_esc else []), library="np")
        decay_add = np.asarray(dd[add_br]) - np.asarray(dd["E_decay_dDep"])
        if has_esc:
            decay_add = decay_add - np.asarray(dd["E_decay_dEsc"])   # remove escaping-daughter energy too
        else:
            print("  (note: no E_decay_dEsc branch — rerun the module so the EndE recovery is neutrino-only;")
            print("   without it, escaping decay daughters are double-counted against the escape term.)")
        lo, hi = np.percentile(decay_add, [1, 99])
        tag = "EndE-dDep-dEsc = neutrino" if (has_esc and not args.decay_rest_mass) else add_br
        print(f"[decay correction ({tag}): mean {decay_add.mean()*1000:+.1f} MeV/event, "
              f"1-99% {lo*1000:+.0f}..{hi*1000:+.0f} MeV]")
        if decay_add.min() < -0.2:
            print(f"  !! WARNING: large NEGATIVE recovery (min {decay_add.min()*1000:.0f} MeV) — "
                  f"expected for the rest-mass variant (undercounts in-flight by the parent KE).")
    edep_all = edep_raw + decay_add       # deposit the CORRECTION builds on (E_dep')
    cc_all   = ev["CCNC"].astype(int)
    cont_all = ev["Contained"].astype(bool) if has_cont else np.ones(len(enu_all), bool)

    # per-event unique key
    ekey_all = np.stack([ev["Run"], ev["SubRun"], ev["Event"]], axis=1)
    ekey_tuple = [tuple(int(x) for x in r) for r in ekey_all]

    sel = (enu_all > 0) & (cc_all == 0)
    if not args.keep_uncontained and has_cont:
        sel &= cont_all
    print(f"[{int(sel.sum())} CC{' contained' if (has_cont and not args.keep_uncontained) else ''} events]")

    # leakage cut — removes genuine high-leakage events the EdgeE containment cut misses
    # (neutral/neutron escape deposits nothing near the edges, so those events pass Contained
    # yet leak most of their energy; they are not reconstructible and poison the channel means).
    if args.max_escapefrac is not None:
        if not has_escf:
            sys.exit("--max-escapefrac needs the EscapeFrac branch on EventTree.")
        escfrac_all = np.asarray(ev["EscapeFrac"], dtype=float)
        before = int(sel.sum())
        leak = sel & (escfrac_all > args.max_escapefrac)
        sel &= escfrac_all <= args.max_escapefrac
        print(f"[leakage cut EscapeFrac <= {args.max_escapefrac:g}: dropped {int(leak.sum())} of {before} "
              f"events ({100*int(leak.sum())/max(before,1):.1f}%); {int(sel.sum())} remain]")

    # build, per selected event, its channel label(s) and the correction target
    if args.level == "event":
        # one channel per event = topology string
        cap = lambda k: min(int(k), 2)
        lepname = {13: "mu", 11: "e"}
        def topo_key(i):
            lp = int(ev["Topo_lepton"][i])
            return (f"{lepname.get(lp,'nolep')}."
                    f"{cap(ev['Topo_nP'][i])}p.{cap(ev['Topo_nPi'][i])}pi.{cap(ev['Topo_nShower'][i])}g")
        for b in ["Topo_lepton", "Topo_nP", "Topo_nPi", "Topo_nShower"]:
            if b not in ev:
                sys.exit(f"EventTree missing '{b}' — needed for --level event")
        idx = np.where(sel)[0]
        chan_of_event = {i: [topo_key(i)] for i in idx}      # event -> [one channel]

    else:  # vertex
        vm = find_tree(f, "VertexMissTree")
        if vm is None:
            sys.exit("no VertexMissTree")
        for b in ["Run", "SubRun", "Event", "Key", "E_missing"]:
            if b not in vm.keys():
                sys.exit(f"VertexMissTree missing '{b}'. If it's Run/SubRun, rebuild the module "
                         f"(the join needs them) and rerun, then point this at the new file.")
        has_inclass = "In_class" in vm.keys()
        cons_on = args.neutron_conservation
        has_shower = "OutShowerKE" in vm.keys()
        extra = []
        if cons_on:
            for b in ["In_KE", "E_binding", "OutChargedKE"]:
                if b not in vm.keys():
                    sys.exit(f"--neutron-conservation needs '{b}' on VertexMissTree; rebuild and rerun.")
            extra = ["In_KE", "E_binding", "OutChargedKE"] + (["OutShowerKE"] if has_shower else [])
        read_v = ["Run", "SubRun", "Event", "Key", "E_missing"] \
                 + (["In_class"] if has_inclass else []) + extra
        v = vm.arrays(read_v, library="np")
        vkey = [tuple(int(x) for x in r) for r in np.stack([v["Run"], v["SubRun"], v["Event"]], axis=1)]
        vKeys = v["Key"]; vMiss = np.asarray(v["E_missing"])
        vInStr = None
        if has_inclass:
            vInStr = np.array([s.decode() if isinstance(s, (bytes, bytearray)) else str(s)
                               for s in v["In_class"]])
        # per-vertex conservation prediction (GeV), clamped >= 0
        inKE = None; cons_pred = None
        if cons_on:
            inKE = np.asarray(v["In_KE"], float)
            cp = inKE - np.asarray(v["E_binding"], float) - np.asarray(v["OutChargedKE"], float)
            if has_shower:
                cp = cp - np.asarray(v["OutShowerKE"], float)
            cons_pred = np.maximum(cp, 0.0)
            if not has_shower:
                print("  (note: no OutShowerKE branch — conservation estimate is NOT subtracting EM "
                      "showers, so it will over-predict; rebuild+rerun with OutShowerKE for the clean version)")
        # neutron-incoming vertices are invisible (a neutron leaves no track) -> drop by default
        drop_n = (not args.keep_neutron_incoming)
        if drop_n and not has_inclass:
            print("  (note: no In_class branch — can't drop neutron-incoming vertices; "
                  "rebuild the module or pass --keep-neutron-incoming)")
            drop_n = False
        def _incoming_is_neutron(j):
            if vInStr is not None:
                return vInStr[j] == "n"
            return str(vKeys[j]).startswith("n->")      # fallback: parse the Key
        # group (channel, miss, cons) by event; cons = conservation pred at observable p->/pi-> vertices
        from collections import defaultdict
        vtx_by_event = defaultdict(list)
        n_drop = 0
        for j in range(len(vMiss)):
            if drop_n and _incoming_is_neutron(j):
                n_drop += 1; continue
            cv = np.nan
            if cons_on and vInStr is not None and vInStr[j] in ("p", "pi") and inKE[j] >= 0:
                cv = cons_pred[j]
            vtx_by_event[vkey[j]].append((str(vKeys[j]), vMiss[j], cv))
        if drop_n:
            print(f"[dropped {n_drop} neutron-incoming vertices (invisible n->... ); "
                  f"{len(vMiss)-n_drop} kept]")
        if cons_on:
            n_cons = sum(1 for lst in vtx_by_event.values() for t in lst if np.isfinite(t[2]))
            print(f"[conservation estimate applied at {n_cons} observable p->/pi-> vertices; "
                  f"nu-> and the rest keep the channel mean]")
        sel_set = set(ekey_tuple[i] for i in np.where(sel)[0])
        idx = np.array([i for i in np.where(sel)[0] if ekey_tuple[i] in vtx_by_event])
        chan_of_event = {i: [t[0] for t in vtx_by_event[ekey_tuple[i]]] for i in idx}
        miss_of_event_vtx = {i: [t[1] for t in vtx_by_event[ekey_tuple[i]]] for i in idx}
        cons_of_event = {i: [t[2] for t in vtx_by_event[ekey_tuple[i]]] for i in idx}

    if len(idx) < 20:
        sys.exit(f"only {len(idx)} usable events — too few.")

    enu  = enu_all[idx]
    edep = edep_all[idx]        # E_dep' (raw + decay recovery) — what the correction builds on
    edepr = edep_raw[idx]       # raw deposit — the deposited-only baseline
    pos  = {ev_i: j for j, ev_i in enumerate(idx)}          # event index -> row in enu/edep

    # ---------- build the channel lookup (mean, std, N) over a set of events ----------
    def build_lookup(train_events):
        """channel -> (mean_miss, std_miss, N, mean_r, std_r), r = E_dep/E_nu (event level).
        Also returns global fallbacks (gmiss, gr) for channels too rare to trust."""
        from collections import defaultdict
        vals = defaultdict(list); rvals = defaultdict(list)
        if args.level == "event":
            for i in train_events:
                c = chan_of_event[i][0]
                vals[c].append(enu_all[i] - edep_all[i])
                if enu_all[i] > 0:
                    rvals[c].append(edep_all[i] / enu_all[i])
        else:
            for i in train_events:
                for c, m in zip(chan_of_event[i], miss_of_event_vtx[i]):
                    vals[c].append(m)
        look = {}
        for c, arr in vals.items():
            a = np.asarray(arr)
            r = np.asarray(rvals[c]) if rvals.get(c) else np.array([np.nan])
            look[c] = (a.mean(), a.std(), len(a), np.nanmean(r), np.nanstd(r))
        gmiss = np.mean([enu_all[i] - edep_all[i] for i in train_events]) if args.level == "event" \
                else np.mean([m for i in train_events for m in miss_of_event_vtx[i]])
        gr = np.mean([edep_all[i] / enu_all[i] for i in train_events if enu_all[i] > 0])
        return look, gmiss, gr

    def correct(events, look, gmiss, gr):
        """E_reco and predicted FRACTIONAL sigma per event, for `events`, per --mode."""
        reco = np.empty(len(events)); sigf = np.empty(len(events))
        for n, i in enumerate(events):
            if args.mode == "additive":
                corr = 0.0; var = 0.0
                cons_list = cons_of_event[i] if args.neutron_conservation else None
                for vj, c in enumerate(chan_of_event[i]):
                    if cons_list is not None and np.isfinite(cons_list[vj]):
                        corr += cons_list[vj]          # conservation estimate at an observable vertex
                        continue
                    m, s, N, mr, sr = look.get(c, (gmiss, 0.0, 0, gr, 0.0))
                    if N < args.min_count:
                        m = gmiss
                    corr += m; var += s * s
                reco[n] = edep_all[i] + corr
                sigf[n] = (np.sqrt(var) / enu_all[i]) if enu_all[i] > 0 else 0.0
            else:  # multiplicative: scale E_dep by the channel's 1/<E_dep/E_nu> (event level)
                c = chan_of_event[i][0]
                m, s, N, mr, sr = look.get(c, (gmiss, 0.0, 0, gr, 0.0))
                if N < args.min_count or not np.isfinite(mr) or mr <= 1e-6:
                    mr, sr = gr, sr
                reco[n] = edep_all[i] / mr
                sigf[n] = (sr / mr) if mr > 1e-6 else 0.0
        return reco, sigf

    # ---------- in-sample (what was asked) ----------
    look_all, gmiss_all, gr_all = build_lookup(idx)
    reco_in, sigf_in = correct(idx, look_all, gmiss_all, gr_all)
    fr_dep  = (edepr - enu) / enu     # deposited-only baseline uses the RAW deposit, never E_dep'
    fr_in   = (reco_in - enu) / enu

    # ---------- k-fold cross-validated (the honest number) ----------
    folds = kfold_indices(len(idx), args.kfolds)
    reco_cv = np.empty(len(idx)); sigf_cv = np.empty(len(idx))
    for fold in folds:
        test_events = idx[fold]
        fold_set = set(int(x) for x in fold)
        train_events = np.array([e for e in idx if pos[e] not in fold_set])
        look, gmiss, gr = build_lookup(train_events)
        r, s = correct(test_events, look, gmiss, gr)
        for n, e in enumerate(test_events):
            reco_cv[pos[e]] = r[n]; sigf_cv[pos[e]] = s[n]
    fr_cv = (reco_cv - enu) / enu

    # ---------- report ----------
    print(f"\n============  NEUTRINO ENERGY RECONSTRUCTION  ({args.level}-level, {args.mode})  ============")
    print(f"channels: {len(look_all)}    events: {len(idx)}    <E_nu> = {enu.mean():.2f} GeV")
    print(f"\n  {'':28s}{'bias':>8}{'RMS':>10}{'68%':>9}")
    print(f"  deposited only (no corr)   : {stat_line(fr_dep)}")
    print(f"  + channel {args.mode[:4]}, in-sample : {stat_line(fr_in)}")
    print(f"  + channel {args.mode[:4]}, {args.kfolds}-fold CV : {stat_line(fr_cv)}   <-- the honest number")

    # does the per-channel spread predict the achieved resolution? (fractional)
    print(f"\n  achieved resolution (CV)   : {fr_cv.std()*100:6.2f}%")
    print(f"  predicted from channel std : {sigf_cv.mean()*100:6.2f}%  (mean per-channel spread)")
    print(f"    -> close agreement means the channel spread explains the resolution; achieved")
    print(f"       >> predicted means extra scatter beyond the channels (closure residual, or the")
    print(f"       correction's functional form is wrong inside the channel, e.g. energy scaling).")

    if args.csv:
        import csv as _csv
        rows = sorted(look_all.items(), key=lambda kv: -kv[1][2])
        with open(args.csv, "w", newline="") as fh:
            w = _csv.writer(fh)
            w.writerow(["channel", "N", "mean_miss_MeV", "std_miss_MeV", "mean_Edep/Enu", "k_mult"])
            for c, (m, s, N, mr, sr) in rows:
                kmul = (1.0 / mr) if (np.isfinite(mr) and mr > 1e-6) else float("nan")
                w.writerow([c, N, f"{m*1000:.2f}", f"{s*1000:.2f}",
                            f"{mr:.4f}" if np.isfinite(mr) else "", f"{kmul:.4f}" if np.isfinite(kmul) else ""])
        print(f"\nwrote correction table -> {args.csv}")

    if args.event_csv:
        import csv as _csv
        with open(args.event_csv, "w", newline="") as fh:
            w = _csv.writer(fh)
            w.writerow(["Run", "SubRun", "Event", "N_vtx",
                        "E_true_GeV", "E_dep_GeV", "E_reco_GeV", "frac_dep", "frac_reco"])
            for j in range(len(idx)):
                orig = int(idx[j])
                rr, ss, ee = ekey_tuple[orig]
                nvtx = len(chan_of_event.get(orig, []))   # # vertices (vertex mode) or 1 (event mode)
                w.writerow([rr, ss, ee, nvtx,
                            f"{enu[j]:.5f}", f"{edepr[j]:.5f}", f"{reco_cv[j]:.5f}",
                            f"{fr_dep[j]:.5f}", f"{fr_cv[j]:.5f}"])   # E_dep raw; frac_reco is the CV value
        print(f"wrote per-event reco -> {args.event_csv}")

    # ---------- plots ----------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("(matplotlib missing — no plot)"); return

    fig, ax = plt.subplots(2, 2, figsize=(13, 10))

    a = ax[0, 0]
    bins = np.linspace(-0.5, 0.3, 60)
    a.hist(fr_dep, bins=bins, color="#C44E52", alpha=0.55,
           label=f"deposited only (RMS {fr_dep.std()*100:.1f}%)")
    a.hist(fr_cv, bins=bins, color="#2C5FAA", alpha=0.7,
           label=f"channel-corrected CV (RMS {fr_cv.std()*100:.1f}%)")
    a.axvline(0, ls="--", color="0.4")
    a.set_xlabel("(E_reco - E_true)/E_true"); a.set_ylabel("events")
    a.set_title("(a) reconstruction resolution"); a.legend()

    a = ax[0, 1]
    hi = np.percentile(enu, 99)
    a.hist2d(enu, reco_cv, bins=60, range=[[0, hi], [0, hi]], cmap="Blues")
    a.plot([0, hi], [0, hi], "r--", lw=1)
    a.set_xlabel("true E_nu (GeV)"); a.set_ylabel("reconstructed E_nu (GeV)")
    a.set_title("(b) reco vs true (channel avg, CV)")

    # resolution & bias vs E_nu
    nb = 10
    be = np.linspace(enu.min(), np.percentile(enu, 97), nb + 1)
    bc, rms_dep, rms_cv, bias_cv = [], [], [], []
    for i in range(nb):
        m = (enu >= be[i]) & (enu < be[i+1])
        if m.sum() >= 15:
            bc.append(0.5*(be[i]+be[i+1]))
            rms_dep.append(fr_dep[m].std()*100); rms_cv.append(fr_cv[m].std()*100)
            bias_cv.append(fr_cv[m].mean()*100)
    a = ax[1, 0]
    a.plot(bc, rms_dep, "o-", color="#C44E52", label="deposited only")
    a.plot(bc, rms_cv, "o-", color="#2C5FAA", label="channel-corrected")
    a.set_xlabel("true E_nu (GeV)"); a.set_ylabel("resolution RMS (%)")
    a.set_title("(c) resolution vs energy"); a.legend()

    a = ax[1, 1]
    a.axhline(0, ls="--", color="0.4")
    a.plot(bc, bias_cv, "o-", color="#2C5FAA")
    a.set_xlabel("true E_nu (GeV)"); a.set_ylabel("bias (%)")
    a.set_title("(d) residual bias vs energy (channel-corrected)")

    plt.tight_layout(); plt.savefig(args.plot, dpi=120)
    print(f"wrote {args.plot}")


if __name__ == "__main__":
    main()
