#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Plot colors of 3I on color-color diagram.
"""
import os 
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt  

from I_common import (
    mycolor, mygrid, I_mags_Bolin, I_colors_Bolin, I_mags_Kareta, I_colors_Kareta,
    stype2colmark, mymathtext, SDSS2PS_mag, SDSS2PS_col, adderr_series, adderr)
from hoya.core import extract_hoya, DBPATH


def rad_vs_col(df, rad_nominal, out, width):
    """Plot radius vs. color.

    Parameters
    ----------
    df : pandas.DataFrame
        input dataframe
    rad_nominal : int
        nominal radius in pix
    out : str
        output filename
    width : float
        y width of the plot
    """

    df_nom = df[df["radius"] == rad_nominal] 
    df = df[df["radius"] != rad_nominal] 

    fig = plt.figure(figsize=(16, 6))
    ax1 = fig.add_axes([0.1, 0.15, 0.25, 0.8])
    ax2 = fig.add_axes([0.4, 0.15, 0.25, 0.8])
    ax3 = fig.add_axes([0.7, 0.15, 0.25, 0.8])
    mygrid(ax1)
    mygrid(ax2)
    mygrid(ax3)

    ax1.set_xlabel("Aperture radius [arcsec]")
    ax2.set_xlabel("Aperture radius [arcsec]")
    ax3.set_xlabel("Aperture radius [arcsec]")
    ax1.set_ylabel("Color [mag]")
    ax1.errorbar(
        df["radius_arcsec"], df["g_r"], df["g_rerr"], fmt="o", color="black", label="g-r")
    ax2.errorbar(
        df["radius_arcsec"], df["r_i"], df["r_ierr"], fmt="o", color="black", label="r-i")
    ax3.errorbar(
        df["radius_arcsec"], df["r_z"], df["r_zerr"], fmt="o", color="black", label="r-z")

    ax1.errorbar(
        df_nom["radius_arcsec"], df_nom["g_r"], df_nom["g_rerr"], fmt="*", 
        ms=15, color=mycolor[0], label="Nominal g-r")
    ax2.errorbar(
        df_nom["radius_arcsec"], df_nom["r_i"], df_nom["r_ierr"], fmt="*", 
        ms=15, color=mycolor[0], label="Nominal r-i")
    ax3.errorbar(
        df_nom["radius_arcsec"], df_nom["r_z"], df_nom["r_zerr"], fmt="*", 
        ms=15, color=mycolor[0], label="Nominal r-z")

    for ax in [ax1, ax2, ax3]:
        ymin, ymax = ax.get_ylim()
        ymean = (ymin + ymax)/2
        ymin1 = ymean + 0.5*width
        ymax1 = ymean - 0.5*width
        ax.set_ylim([ymin1, ymax1])
    
    loc = "upper left"
    fs = 20
    ax1.legend(loc=loc, fontsize=fs)
    ax2.legend(loc=loc, fontsize=fs)
    ax3.legend(loc=loc, fontsize=fs)
    fig.savefig(out)
    plt.close()


if __name__ == "__main__":
    parser = ap(description="Plot color-color diagram.")
    parser.add_argument(
        "col_I", type=str, 
        help="Colors of 3I")
    parser.add_argument(
        "col_P", type=str, 
        help="Colors of Popescu")
    parser.add_argument(
        "radius", type=int, 
        help="Aperture radius in pix")
    parser.add_argument(
        "--Popescu", action="store_true", default=None,
        help="Plot colors of Popescu")
    parser.add_argument(
        "--magtype", type=str, default="PS",
        help="Magnitude system")
    parser.add_argument(
        "--Bolin", action="store_true", default=None,
        help="Plot colors from Bolin+2025")
    parser.add_argument(
        "--Kareta", action="store_true", default=None,
        help="Plot colors from Kareta+2025")
    parser.add_argument(
        "--outtype", type=str, default="pdf",
        help="format of output figures")
    parser.add_argument(
        "--outdir", type=str, default="fig",
        help="output directory")
    args = parser.parse_args()
    
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)
   
    # arcsec/pix
    p_scale = 0.35
    col = ["obj", "stype", "H"]
    magtype = args.magtype
    stype_use = ["S", "V", "X", "K", "L", "C", "B", "D", "A"]
    col_Seimei, s_Seimei, mark_Seimei = mycolor[1], 700, "*"
    rad_pix = args.radius
    rad_arcsec = args.radius*p_scale
    print(f"Use data with aperture radius of {rad_pix} pix = {rad_arcsec:.3f} arcsec")

    df_I = pd.read_csv(args.col_I, sep=" ")
    df_P = pd.read_csv(args.col_P, sep=" ")
    df_I["radius_arcsec"] = df_I["radius"]* p_scale
    df_P["radius_arcsec"] = df_P["radius"]* p_scale

    # Radius vs. color ========================================================
    # 3I
    out = f"I_fig_radcol.{args.outtype}"
    out = os.path.join(args.outdir, out)
    rad_vs_col(df_I, rad_pix, out, width=0.3)
    # Popescu
    out = f"P_fig_radcol.{args.outtype}"
    out = os.path.join(args.outdir, out)
    rad_vs_col(df_I, rad_pix, out, width=0.3)
    # Radius vs. color ========================================================


    # Read Seimei colors ======================================================
    df_I_r = df_I[df_I["radius"] == args.radius]
    g_r_I, g_rerr_I = df_I_r["g_r"].iloc[0], df_I_r["g_rerr"].iloc[0]
    r_i_I, r_ierr_I = df_I_r["r_i"].iloc[0], df_I_r["r_ierr"].iloc[0]
    r_z_I, r_zerr_I = df_I_r["r_z"].iloc[0], df_I_r["r_zerr"].iloc[0]
    # Calculate i_z
    i_z_I    = r_z_I - r_i_I
    i_zerr_I = adderr(r_zerr_I, r_ierr_I)

    print(f" 3I/ATLAS Seimei Pan-STARRS (original)")
    print(f"    g-r {g_r_I:.3f}+-{g_rerr_I:.3f}")
    print(f"    r-i {r_i_I:.3f}+-{r_ierr_I:.3f}")
    print(f"    i-z {i_z_I:.3f}+-{i_zerr_I:.3f}")
    print(f"    r-z {r_z_I:.3f}+-{r_zerr_I:.3f}")
    print("")

    df_P_r = df_P[df_P["radius"] == args.radius]
    g_r_P, g_rerr_P = df_P_r["g_r"].iloc[0], df_P_r["g_rerr"].iloc[0]
    r_i_P, r_ierr_P = df_P_r["r_i"].iloc[0], df_P_r["r_ierr"].iloc[0]
    r_z_P, r_zerr_P = df_P_r["r_z"].iloc[0], df_P_r["r_zerr"].iloc[0]

    print(f" Popescu Seimei Pan-STARRS (original)")
    print(f"    g-r {g_r_P:.3f}+-{g_rerr_P:.3f}")
    print(f"    r-i {r_i_P:.3f}+-{r_ierr_P:.3f}")
    print(f"    r-z {r_z_P:.3f}+-{r_zerr_P:.3f}")
    print("")
    # Read Seimei colors ======================================================


    # Other papers ============================================================
    # 1. Bolin+2025, SDSS
    g, gerr = I_mags_Bolin["g"], I_mags_Bolin["gerr"]
    r, rerr = I_mags_Bolin["r"], I_mags_Bolin["rerr"]
    i, ierr = I_mags_Bolin["i"], I_mags_Bolin["ierr"]
    z, zerr = I_mags_Bolin["z"], I_mags_Bolin["zerr"]

    df_B25 = pd.DataFrame(dict(
        g=[g], gerr=[gerr],
        r=[r], rerr=[rerr],
        i=[i], ierr=[ierr],
        z=[z], zerr=[zerr]
        ))
    df_B25 = SDSS2PS_mag(df_B25)
    g_r_B    = (df_B25["g_PS"] - df_B25["r_PS"])[0]
    g_rerr_B = (adderr_series(df_B25["gerr_PS"], df_B25["rerr_PS"]))[0]
    r_i_B    = (df_B25["r_PS"] - df_B25["i"])[0]
    r_ierr_B = (adderr_series(df_B25["rerr_PS"], df_B25["ierr_PS"]))[0]
    r_z_B    = (df_B25["r_PS"] - df_B25["z_PS"])[0]
    r_zerr_B = (adderr_series(df_B25["rerr_PS"], df_B25["zerr_PS"]))[0]

    print(f" 3I/ATLAS Bolin Pan-STARRS (converted, from mag)")
    print(f"    g-r {g_r_B:.3f}+-{g_rerr_B:.3f}")
    print(f"    r-i {r_i_B:.3f}+-{r_ierr_B:.3f}")
    print(f"    r-z {r_z_B:.3f}+-{r_zerr_B:.3f}")
    print("")
    

    g_r, g_rerr = I_colors_Bolin["g_r"], I_colors_Bolin["g_rerr"]
    r_i, r_ierr = I_colors_Bolin["r_i"], I_colors_Bolin["r_ierr"]
    i_z, i_zerr = I_colors_Bolin["i_z"], I_colors_Bolin["i_zerr"]
    df_B25_col_SDSS = pd.DataFrame(dict(
        g_r_SDSS=[g_r], g_rerr_SDSS=[g_rerr],
        r_i_SDSS=[r_i], r_ierr_SDSS=[r_ierr],
        i_z_SDSS=[i_z], i_zerr_SDSS=[i_zerr],
        ))
    df_B25_col_PS = SDSS2PS_col(df_B25_col_SDSS)
    g_r_B    = df_B25_col_PS["g_r_PS"].iloc[0]
    g_rerr_B = df_B25_col_PS["g_rerr_PS"].iloc[0]
    r_i_B    = df_B25_col_PS["r_i_PS"].iloc[0]
    r_ierr_B = df_B25_col_PS["r_ierr_PS"].iloc[0]
    i_z_B    = df_B25_col_PS["i_z_PS"].iloc[0]
    i_zerr_B = df_B25_col_PS["i_zerr_PS"].iloc[0]
    # Calculate r_z
    r_z_B    = r_i_B + i_z_B
    r_zerr_B = adderr(r_ierr_B, i_zerr_B)
    print(f" 3I/ATLAS Bolin Pan-STARRS (converted, from colors)")
    print(f"    g-r {g_r_B:.3f}+-{g_rerr_B:.3f}")
    print(f"    r-i {r_i_B:.3f}+-{r_ierr_B:.3f}")
    print(f"    i-z {i_z_B:.3f}+-{i_zerr_B:.3f}")
    print(f"    r-z {r_z_B:.3f}+-{r_zerr_B:.3f}")
    print("")

    # 2. Kareta+2025, Kareta
    print(f" 3I/ATLAS Kareta Pan-STARRS (original)")
    g_i_K, g_ierr_K = I_colors_Kareta["g_i"], I_colors_Kareta["g_ierr"]
    print(f"    g-i {g_i_K:.3f}+-{g_ierr_K:.3f}")
    print("")
    
    # ZTF colors in Seligman+2025
    # Maybe Pan-STARRS (Masci+2018)
    # TODO: Check 
    # 2025-05-22
    g_r_Z1 = 0.42
    # 2025-06-18
    g_r_Z2 = 0.44
    # Other papers ============================================================


    # Read SDSS MOC data ======================================================
    df_S21 = extract_hoya(DBPATH, "Sergeyev2021")
    # Select high possibility objects
    p_th = 0.8
    df_S21 = df_S21[df_S21["p_stype"].astype(float) > p_th]
    print(f"N_SDSS = {len(df_S21)} (p > {p_th})")

    if magtype=="PS":
        # Change system to PS system
        print("Photometric results are in the Pan-STARRS.")
        print("Convert colors in the SDSS to those in PS.")
        df_S21 = SDSS2PS_mag(df_S21)

        # Rename
        df_S21 = df_S21.rename(
          columns={"g":"g_temp", "gerr":"gerr_temp", 
                   "r":"r_temp", "rerr":"rerr_temp", 
                   "i":"i_temp", "ierr":"ierr_temp", 
                   "z":"z_temp", "zerr":"zerr_temp"})
        df_S21 = df_S21.rename(
          columns={"g_PS":"g", "gerr_PS":"gerr", 
                   "r_PS":"r", "rerr_PS":"rerr", 
                   "i_PS":"i", "ierr_PS":"ierr", 
                   "z_PS":"z", "zerr_PS":"zerr"})
    
    # Popescu
    df_P_SDSS = df_S21[df_S21["obj"]=="Popescu"]
    g_r_P_SDSS = df_P_SDSS["g"] - df_P_SDSS["r"]
    g_rerr_P_SDSS = adderr(df_P_SDSS["gerr"], df_P_SDSS["rerr"])
    r_i_P_SDSS = df_P_SDSS["r"] - df_P_SDSS["i"]
    r_ierr_P_SDSS = adderr(df_P_SDSS["rerr"], df_P_SDSS["ierr"])
    r_z_P_SDSS = df_P_SDSS["r"] - df_P_SDSS["z"]
    r_zerr_P_SDSS = adderr(df_P_SDSS["rerr"], df_P_SDSS["zerr"])
    print(f" Popescu SDSS in the  Pan-STARRS system (Converted)")
    print(f"    g-r {g_r_P_SDSS.values[0]:.3f}+-{g_rerr_P_SDSS:.3f}")
    print(f"    r-i {r_i_P_SDSS.values[0]:.3f}+-{r_ierr_P_SDSS:.3f}")
    print(f"    r-z {r_z_P_SDSS.values[0]:.3f}+-{r_zerr_P_SDSS:.3f}")
    print("")
    # Read SDSS MOC data ======================================================


    # 1. g-r vs. r-i/r-z ======================================================
    mgtp = "Pan-STARRS"
    fig = plt.figure(figsize=(8, 16))
    ax_u = fig.add_axes([0.15, 0.56, 0.80, 0.40])
    ax_l = fig.add_axes([0.15, 0.10, 0.80, 0.40])
    ax_u.set_ylabel(f"$r-i$  ({mgtp})")
    ax_l.set_xlabel(f"$g-r$  ({mgtp})")
    ax_l.set_ylabel(f"$r-z$  ({mgtp})")

    # Set filters width = 0.8
    c1min, c1max = 0.2, 1.0
    c2min, c2max = -0.2, 0.6
    c3min, c3max = -0.3, 0.5
    print(f"Color range")
    print(f"    c1 = {c1min}--{c1max}")
    print(f"    c2 = {c2min}--{c2max}")
    print(f"    c3 = {c3min}--{c3max}")

    ax_u.set_xlim([c1min, c1max])
    ax_u.set_ylim([c2min, c2max])
    ax_l.set_xlim([c1min, c1max])
    ax_l.set_ylim([c3min, c3max])
    
    mymathtext()

    # Plot objects in SDSSMOC
    for idx, stype in enumerate(stype_use):
        df_type = df_S21[df_S21["stype"]==stype]
        df_type = df_type.reset_index(drop=True)
        print(f"  {stype}-complex N={len(df_type)}")
        if len(df_type) == 0:
            continue
        col, mark = stype2colmark(stype)
    
        s_SDSS, lw_SDSS = 50, 1
        ax_u.scatter(
            df_type["g"]-df_type["r"],
            df_type["r"]-df_type["i"],
            color=col, s=s_SDSS, lw=lw_SDSS, marker=mark, facecolor="None",
            edgecolor=col, zorder=-1, label=f"{stype}-complex")
        ax_l.scatter(
            df_type["g"]-df_type["r"],
            df_type["r"]-df_type["z"],
            color=col, s=s_SDSS, lw=lw_SDSS, marker=mark, facecolor="None",
            edgecolor=col, zorder=-1, label=f"{stype}-complex")
    
    # 3I Seimei ===============================================================
    col, mark, ms = mycolor[0], "*", 25
    ax_u.errorbar(
        g_r_I, r_i_I, xerr=g_rerr_I, yerr=r_ierr_I,
        fmt=mark, markerfacecolor=col, 
        markeredgecolor="black", ecolor='black', 
        ms=ms, lw=2, zorder=201, label="3I/ATLAS\n  (This study)",
        )
    ax_l.errorbar(
        g_r_I, r_z_I, xerr=g_rerr_I, yerr=r_zerr_I,
        fmt=mark, markerfacecolor=col, 
        markeredgecolor="black", ecolor='black', 
        ms=ms, lw=2, zorder=201, label="3I/ATLAS\n  (This study)",
        )
    # 3I Seimei ===============================================================


    # 3I Bolin ================================================================
    col, mark, ms = mycolor[1], "o", 20
    ax_u.errorbar(
        g_r_B, r_i_B, xerr=g_rerr_B, yerr=r_ierr_B,
        fmt=mark, markerfacecolor=col, 
        markeredgecolor="black", ecolor='black', 
        ms=ms, lw=2, zorder=201, label="3I/ATLAS\n  (Bolin+2025)",
        )
    ax_l.errorbar(
        g_r_B, r_z_B, xerr=g_rerr_B, yerr=r_zerr_B,
        fmt=mark, markerfacecolor=col, 
        markeredgecolor="black", ecolor='black', 
        ms=ms, lw=2, zorder=201, label="3I/ATLAS\n  (Bolin+2025)",
        )
    # 3I Bolin ================================================================


    # 3I Kareta ===============================================================
    #col, mark = "orange", "^"
    #ax_u.errorbar(
    #    g_r_K, r_i_K, xerr=g_rerr_K, yerr=r_ierr_K,
    #    fmt=mark, markerfacecolor=col, 
    #    markeredgecolor="black", ecolor='black', 
    #    ms=25, lw=2, zorder=201, label="3I/ATLAS\n  (Kareta+2025)",
    #    )
    # No z-band measurements
    # 3I Kareta ===============================================================


    # Popescu colors ==========================================================
    if args.Popescu:
        # Seimei
        col, mark = mycolor[0], "h"
        ax_u.errorbar(
            g_r_P, r_i_P, xerr=g_rerr_P, yerr=r_ierr_P,
            fmt=mark, markerfacecolor=col, 
            markeredgecolor="black", ecolor='black', 
            ms=25, lw=2, zorder=201, label="Popescu\n  (This study)",
            )
        ax_l.errorbar(
            g_r_P, r_z_P, xerr=g_rerr_P, yerr=r_zerr_P,
            fmt=mark, markerfacecolor=col, 
            markeredgecolor="black", ecolor='black', 
            ms=25, lw=2, zorder=201, label="Popescu\n  (This study)",
            )
        # SDSS ================================================================

        col, mark = mycolor[0], "H"
        ax_u.errorbar(
            g_r_P_SDSS, r_i_P_SDSS, xerr=g_rerr_P_SDSS, yerr=r_ierr_P_SDSS,
            fmt=mark, markerfacecolor=col, 
            markeredgecolor="black", ecolor='black', 
            ms=25, lw=2, zorder=201, label="Popescu$\n  (SDSS)",
            )
        ax_l.errorbar(
            g_r_P_SDSS, r_z_P_SDSS, xerr=g_rerr_P_SDSS, yerr=r_zerr_P_SDSS,
            fmt=mark, markerfacecolor=col, 
            markeredgecolor="black", ecolor='black', 
            ms=25, lw=2, zorder=201, label="Popescu\n  (SDSS)",
            )
    # Popescu colors ==========================================================

    for ax in [ax_l, ax_u]:
        ax.legend(loc="upper right", fontsize=14).get_frame().set_alpha(1.0)

    ax_u.text(-0.12, 1.05, "(a)", size=22, transform=ax_u.transAxes)
    ax_l.text(-0.12, 1.05, "(b)", size=22, transform=ax_l.transAxes)
    
    out = f"I_fig_cc_{magtype}.{args.outtype}"
    out = os.path.join(args.outdir, out)
    fig.savefig(out)
    plt.close()
    # 1. g-r vs. r-i/r-z ======================================================
