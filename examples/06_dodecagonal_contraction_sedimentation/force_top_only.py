#!/usr/bin/env python3

# python3 force_top_only.py log.pimpleLYJHFDIBFoam --top-n 100 > force_summary.txt

import re
import sys
import argparse
import heapq
from collections import defaultdict

num = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?"

def make_pat(name):
    return re.compile(r"\b" + re.escape(name) + r"\s+(" + num + r")")

body_pat = re.compile(r"\bbody\s+(\d+)")
patch_pat = re.compile(r"\bpatches\s+(.*?)\s+contactVolume\s+")

pats = {
    "contactVolumeOverH3": make_pat("contactVolumeOverH3"),
    "contactAreaOverH2": make_pat("contactAreaOverH2"),
    "magFNe": make_pat("magFNe"),
    "magFNd": make_pat("magFNd"),
    "magFt": make_pat("magFt"),
    "magFA": make_pat("magFA"),
    "magFtotal": make_pat("magFtotal"),
    "fPerVol": make_pat("fPerVol"),
    "Lc": make_pat("Lc"),
    "Vn": make_pat("Vn"),
}

flag_pats = {
    "flag_bigVol": re.compile(r"\bflag_bigVol\s+([01])"),
    "flag_bigArea": re.compile(r"\bflag_bigArea\s+([01])"),
    "flag_smallLc": re.compile(r"\bflag_smallLc\s+([01])"),
    "flag_bigFPerVol": re.compile(r"\bflag_bigFPerVol\s+([01])"),
}

def get_float(line, key, default=0.0):
    m = pats[key].search(line)
    if not m:
        return default
    return float(m.group(1))

def get_flag(line, key):
    m = flag_pats[key].search(line)
    if not m:
        return 0
    return int(m.group(1))

def push_top(heap, top_n, score, rec, serial):
    if top_n <= 0:
        return

    item = (score, serial, rec)

    if len(heap) < top_n:
        heapq.heappush(heap, item)
    elif score > heap[0][0]:
        heapq.heapreplace(heap, item)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("logfile")
    ap.add_argument("--top-n", type=int, default=100)
    ap.add_argument("--min-force", type=float, default=0.0)
    ap.add_argument("--min-fpervol", type=float, default=0.0)
    ap.add_argument("--min-vol-h3", type=float, default=0.0)
    args = ap.parse_args()

    body_sum = defaultdict(lambda: {
        "n": 0,
        "maxF": 0.0,
        "maxFNe": 0.0,
        "maxFNd": 0.0,
        "maxFt": 0.0,
        "maxFA": 0.0,
        "maxFPerVol": 0.0,
        "maxVolH3": 0.0,
        "maxAreaH2": 0.0,
        "minLc": None,
        "bigVol": 0,
        "bigArea": 0,
        "smallLc": 0,
        "bigFPerVol": 0,
        "domFNe": 0,
        "domFNd": 0,
        "domFt": 0,
        "domFA": 0,
        "lastPatch": "",
        "lastLine": 0,
    })

    top_force = []
    top_fpervol = []
    top_vol = []

    n_lines = 0
    n_force_lines = 0
    n_selected = 0
    n_parse_fail = 0
    serial = 0

    with open(args.logfile, "r", errors="ignore") as f:
        for lineno, line in enumerate(f, 1):
            n_lines = lineno

            if "WALLDBG_FORCE_EXPLODE" not in line:
                continue

            n_force_lines += 1

            bm = body_pat.search(line)
            if not bm:
                n_parse_fail += 1
                continue

            body = int(bm.group(1))

            pm = patch_pat.search(line)
            patch = pm.group(1).strip() if pm else "NA"

            magF = get_float(line, "magFtotal")
            fPerVol = get_float(line, "fPerVol")
            volH3 = get_float(line, "contactVolumeOverH3")
            areaH2 = get_float(line, "contactAreaOverH2")
            lc = get_float(line, "Lc", -1.0)

            if magF < args.min_force:
                continue
            if fPerVol < args.min_fpervol:
                continue
            if volH3 < args.min_vol_h3:
                continue

            n_selected += 1

            magFNe = get_float(line, "magFNe")
            magFNd = get_float(line, "magFNd")
            magFt = get_float(line, "magFt")
            magFA = get_float(line, "magFA")

            comps = [
                ("FNe", magFNe),
                ("FNd", magFNd),
                ("Ft", magFt),
                ("FA", magFA),
            ]
            dom, dom_val = max(comps, key=lambda x: x[1])

            bigVol = get_flag(line, "flag_bigVol")
            bigArea = get_flag(line, "flag_bigArea")
            smallLc = get_flag(line, "flag_smallLc")
            bigFPerVol = get_flag(line, "flag_bigFPerVol")

            c = body_sum[body]
            c["n"] += 1
            c["maxF"] = max(c["maxF"], magF)
            c["maxFNe"] = max(c["maxFNe"], magFNe)
            c["maxFNd"] = max(c["maxFNd"], magFNd)
            c["maxFt"] = max(c["maxFt"], magFt)
            c["maxFA"] = max(c["maxFA"], magFA)
            c["maxFPerVol"] = max(c["maxFPerVol"], fPerVol)
            c["maxVolH3"] = max(c["maxVolH3"], volH3)
            c["maxAreaH2"] = max(c["maxAreaH2"], areaH2)

            if lc >= 0.0 and (c["minLc"] is None or lc < c["minLc"]):
                c["minLc"] = lc

            c["bigVol"] += bigVol
            c["bigArea"] += bigArea
            c["smallLc"] += smallLc
            c["bigFPerVol"] += bigFPerVol
            c["lastPatch"] = patch
            c["lastLine"] = lineno

            if dom == "FNe":
                c["domFNe"] += 1
            elif dom == "FNd":
                c["domFNd"] += 1
            elif dom == "Ft":
                c["domFt"] += 1
            else:
                c["domFA"] += 1

            rec = (
                lineno,
                body,
                patch,
                magF,
                fPerVol,
                volH3,
                areaH2,
                lc,
                dom,
                dom_val,
                magFNe,
                magFNd,
                magFt,
                magFA,
                bigVol,
                bigArea,
                smallLc,
                bigFPerVol,
            )

            serial += 1
            push_top(top_force, args.top_n, magF, rec, serial)
            push_top(top_fpervol, args.top_n, fPerVol, rec, serial)
            push_top(top_vol, args.top_n, volH3, rec, serial)

    print("=== FORCE SUMMARY ===")
    print(f"lines_read {n_lines}")
    print(f"force_debug_lines {n_force_lines}")
    print(f"selected_lines {n_selected}")
    print(f"parse_fail {n_parse_fail}")
    print(f"unique_bodies {len(body_sum)}")

    print("\n=== BODY TOP BY maxF ===")
    ranked = sorted(
        body_sum.items(),
        key=lambda kv: (
            kv[1]["maxF"],
            kv[1]["maxFPerVol"],
            kv[1]["maxVolH3"],
            kv[1]["n"],
        ),
        reverse=True,
    )

    for body, c in ranked[:args.top_n]:
        minLc = c["minLc"] if c["minLc"] is not None else -1.0
        print(
            f"body {body} "
            f"n {c['n']} "
            f"maxF {c['maxF']:.8g} "
            f"maxFNe {c['maxFNe']:.8g} "
            f"maxFNd {c['maxFNd']:.8g} "
            f"maxFt {c['maxFt']:.8g} "
            f"maxFA {c['maxFA']:.8g} "
            f"maxFPerVol {c['maxFPerVol']:.8g} "
            f"maxVolH3 {c['maxVolH3']:.8g} "
            f"maxAreaH2 {c['maxAreaH2']:.8g} "
            f"minLc {minLc:.8g} "
            f"flags bigVol {c['bigVol']} bigArea {c['bigArea']} "
            f"smallLc {c['smallLc']} bigFPerVol {c['bigFPerVol']} "
            f"dom FNe {c['domFNe']} FNd {c['domFNd']} "
            f"Ft {c['domFt']} FA {c['domFA']} "
            f"lastPatch {c['lastPatch']} "
            f"lastLine {c['lastLine']}"
        )

    def print_heap(title, heap):
        print(f"\n=== {title} ===")
        for _, _, r in sorted(heap, key=lambda x: x[0], reverse=True):
            (
                lineno,
                body,
                patch,
                magF,
                fPerVol,
                volH3,
                areaH2,
                lc,
                dom,
                dom_val,
                magFNe,
                magFNd,
                magFt,
                magFA,
                bigVol,
                bigArea,
                smallLc,
                bigFPerVol,
            ) = r

            print(
                f"L {lineno} "
                f"body {body} "
                f"patch {patch} "
                f"magF {magF:.8g} "
                f"fPerVol {fPerVol:.8g} "
                f"volH3 {volH3:.8g} "
                f"areaH2 {areaH2:.8g} "
                f"Lc {lc:.8g} "
                f"dominant {dom} "
                f"domVal {dom_val:.8g} "
                f"FNe {magFNe:.8g} "
                f"FNd {magFNd:.8g} "
                f"Ft {magFt:.8g} "
                f"FA {magFA:.8g} "
                f"bigVol {bigVol} "
                f"bigArea {bigArea} "
                f"smallLc {smallLc} "
                f"bigFPerVol {bigFPerVol}"
            )

    print_heap("TOP BY magFtotal", top_force)
    print_heap("TOP BY fPerVol", top_fpervol)
    print_heap("TOP BY contactVolumeOverH3", top_vol)

if __name__ == "__main__":
    main()