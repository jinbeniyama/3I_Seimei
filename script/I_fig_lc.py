#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Create lightcurves of 3I/ATLAS.

Input file is 3I_mag_20250715.txt (Seimei obs.)
with "jd", "mag", "magerr", "frameid", and "band".
"""
import os
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from I_common import bandcolor, bandmark, mygrid


if __name__ == "__main__":
    parser = ap(
        description="Create lightcurves of 3I/ATLAS")
    parser.add_argument(
        "res", type=str,
        help="Photometric results of Seimei TriCCS")
    parser.add_argument(
        "--outdir", type=str, default="fig",
        help="output directory")
    parser.add_argument(
        "--outtype", default="pdf",
        help="format of output figure")
    args = parser.parse_args()

    outdir = args.outdir
    if not os.path.isdir(outdir):
      os.makedirs(outdir)
     
    # LT corrected time calculated based on t_jd_g
    bands = ["g", "r", "i", "z"]

    # Read Seimei data as df
    df = pd.read_csv(args.res, sep=" ")

    # Extracd jd in the r band with frameid, and update jds.
    # This is necesarry to link gri or grz observations.
    frameid_list = sorted(list(set(df.frameid)))
    df_r = df[df["band"] == "r"]
    jd_r_list = []
    for frameid in frameid_list:
        df_r_1 = df_r[df_r["frameid"] == frameid]
        jd_r = np.min(df_r_1["jd"])
        df.loc[df["frameid"] == frameid, "jd"] = jd_r

    JD0 = np.min(df["jd"])
    df["t_min"] = (df["jd"] - JD0)*24.*60.
    arc_min = np.max(df["t_min"])
    print(f"  Observation arc: {arc_min:.2f} min")
    print(f"                   {arc_min/60:.2f} hr")


    # Observed magnitude JD vs. mag
    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_axes([0.15, 0.17, 0.80, 0.80])
    ax.set_xlabel("Elapsed time [min]")
    ax.set_ylabel("Observed magnitude [mag]")
    mygrid(ax)
    
    for b in bands:
        df_b = df[df["band"] == b]
        col, mark = bandcolor[b], bandmark[b]
        ax.errorbar(
            df_b["t_min"], df_b["mag"], yerr=df_b["magerr"],
            fmt=mark,
            mfc="None", mec=col, ecolor=col, ms=12,
            capsize=0, zorder=1, lw=1, label=f"{b}")

    ax.set_ylim([19.0, 17.25])
    ax.legend(fontsize=13, loc="lower left", ncol=4).get_frame().set_alpha(1)


    out = f"I_fig_lc.{args.outtype}"
    out = os.path.join(args.outdir, out)
    fig.savefig(out)
    plt.close()
