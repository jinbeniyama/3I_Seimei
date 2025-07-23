#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Plot 3I's reflectance.
"""
import os 
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt  
from matplotlib.ticker import FormatStrFormatter
import classy

from I_common import (
    mycolor, calc_ref_from_color,adderr_series, mygrid, lam_W18)

# Color 
col_stype = dict(A="green", S="green", Q="blue", D="red", Z="blue")


if __name__ == "__main__":
    parser = ap(description="Plot reflectance.")
    parser.add_argument(
        "col_I", type=str,
        help="Colors of 3I")
    parser.add_argument(
        "radius", type=int,
        help="Aperture radius in pix")
    parser.add_argument(
        "--xr", type=float, nargs=2, default=[0.40, 1.0],
        help="x range")
    parser.add_argument(
        "--yr", type=float, nargs=2, default=[0.8, 1.6],
        help="x range")
    parser.add_argument(
        "--outtype", type=str, default="jpg",
        help="format of output figures")
    parser.add_argument(
        "--outdir", type=str, default="fig",
        help="output directory")
    args = parser.parse_args()
    
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    # Common setting in VIS and NIR
    col_I, mark_I = mycolor[0], "*"
    label_I = "3I/ATLAS (This study)"
    st_template = ["S", "D", "Z"]
    fs = 20
    sys = "PS"
    p_scale = 0.35
    # Normalize at g-band (PS) center
    w_norm_vis = lam_W18("g")[0]

    rad_pix = args.radius
    rad_arcsec = args.radius*p_scale
    print(f"Use data with aperture radius of {rad_pix} pix = {rad_arcsec:.3f} arcsec")

    df = pd.read_csv(args.col_I, sep=" ")
    df = df[df["radius"] == rad_pix].copy()
    df["g_g"] = 0
    df["g_gerr"] = 0
    df["r_g"] = -df["g_r"]
    df["r_gerr"] = df["g_rerr"]
    df["i_g"] = - (df["g_r"] + df["r_i"])
    df["i_gerr"] = adderr_series(df["g_rerr"], df["r_ierr"])
    df["z_g"] = - (df["g_r"] + df["r_z"])
    df["z_gerr"] = adderr_series(df["g_rerr"], df["r_zerr"])

    # Seimei 3I
    df_ref = calc_ref_from_color(df, sys, center="g")

    fig = plt.figure(figsize=(8, 8))  
    ax = fig.add_axes([0.15, 0.15, 0.80, 0.80])
    ax.set_xlabel(r"Wavelength [$\mu$m]")
    ax.set_ylabel("Normalized reflectance")
    mygrid(ax)
    
    w_I = [
        lam_W18("g")[0], lam_W18("r")[0], lam_W18("i")[0], lam_W18("z")[0]
        ]
    werr_I = [
        lam_W18("g")[1]/2, lam_W18("r")[1]/2, lam_W18("i")[1]/2, lam_W18("z")[1]/2]
    ref_I = [
        df_ref["ref_g"].iloc[0], df_ref["ref_r"].iloc[0], 
        df_ref["ref_i"].iloc[0], df_ref["ref_z"].iloc[0], 
        ]
    referr_I = [
        df_ref["referr_g"].iloc[0], df_ref["referr_r"].iloc[0], 
        df_ref["referr_i"].iloc[0], df_ref["referr_z"].iloc[0], 
        ]

    ax.errorbar(
        w_I, ref_I, xerr=werr_I, yerr=referr_I,
        fmt=mark_I, markerfacecolor=col_I, markeredgecolor="black",
        ecolor=col_I, ms=25, zorder=3000, label=label_I,
        )

    # Plot template
    temp = classy.taxonomies.mahlke.load_templates()
    for st in st_template:
        sp = temp[st]
        col = col_stype[st]
        label=f"{st}-type"
        # Normalize at g-band center
        sp.normalize(at=w_norm_vis)
        ax.plot(
            sp.wave, sp.refl, lw=1, alpha=0.1, color="gray")
        ax.fill_between(
            sp.wave, sp.refl-sp.refl_err, sp.refl+sp.refl_err, color=col, alpha=0.2, label=label)

    # Change orders of legends
    # 0,1,2,3 -> 3,0,1,2
    hs,ls = ax.get_legend_handles_labels()
    handles = [hs[3], hs[0], hs[1], hs[2]]
    labels = [ls[3], ls[0], ls[1], ls[2]]
    ax.legend(handles, labels, fontsize=fs, loc="upper left")

    ax.set_yticks(np.arange(0, 1.71, 0.10))
    ax.set_xlim([0.40, 1.00])
    ax.set_ylim([0.80, 1.70])
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    out = f"I_fig_ref.{args.outtype}"
    out = os.path.join(args.outdir, out)
    fig.savefig(out)
    plt.close()
