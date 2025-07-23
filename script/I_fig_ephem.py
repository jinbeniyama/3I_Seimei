#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Create ephemeris plot of 3I. (Not used in the paper)
"""
import os
import datetime
from argparse import ArgumentParser as ap
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.dates import date2num
from matplotlib.ticker import MultipleLocator
from astroquery.jplhorizons import Horizons


def timelabel(ax, ttype):
    """
    Parameters
    ----------
    ax : matplotlib.axes
        axis of the plot
    ttype : str
        hour, day, or month

    Return
    ------
    ax : matplotlib.axes
        updated axis of the plot
    """
    if ttype=="hour":
        ax.xaxis.set_major_locator(HourLocator(byhour=range(0, 24, 1)))
        ax.xaxis.set_major_formatter(DateFormatter('%b-%d %H:%M'))
    if ttype=="day":
        ax.xaxis.set_major_locator(DayLocator())
        ax.xaxis.set_major_formatter(DateFormatter('%b-%d'))
    if ttype=="month":
        ax.xaxis.set_major_locator(MonthLocator())
        ax.xaxis.set_major_formatter(DateFormatter('%b-%d'))
    return ax


if __name__ == "__main__":
    parser = ap(
        description="Create ephemeris plot of 3I/ATLAS")
    parser.add_argument(
        "--outdir", type=str, default="fig",
        help="output directory")
    parser.add_argument(
        "--ttype", type=str, default="week", 
        help="Time in labels (hour, day, week, month)")
    parser.add_argument(
        "--outtype", default="pdf",
        help="format of output figure")
    parser.add_argument(
        "--obtain", action="store_true",
        help="Obtain JPL ehpemeris")
    args = parser.parse_args()

    outdir = args.outdir
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    ephem = f"../data/ephem_3I_202500601to20260131_500.txt"

    # Obtain JPL ephem and save it
    if args.obtain:
        date0 = "2025-06-01"
        date1 = "2026-02-28"
        step = "1h"
        obj = Horizons(
            id="3I", location=500, 
            epochs={'start':date0, 'stop':date1, 'step':step})
        eph = obj.ephemerides()
        eph.write(ephem, format="pandas.csv", sep=" ")

    # Read ephemeris data obtained in JPL/HORIZONS
    df = pd.read_csv(ephem, sep=" ")
    
    # For ax 1
    ymin1_l, ymax1_l = 20, 11
    ymin1_r, ymax1_r = 0, 35
    # For ax 
    ymin2_l, ymax2_l = 1, 5
    ymin2_r, ymax2_r = 1, 5

    col_V, col_a = "black", "gray"
    ls_V, ls_a   = "solid", "dashed"

    key_mag = "Tmag"
 
    fig = plt.figure(figsize=(16, 8))
    # V-mag (left) & alpha (right)
    ax1 = fig.add_axes([0.08,  0.60, 0.85, 0.34])
    # delta (left) and r (right)
    ax2 = fig.add_axes([0.08,  0.18, 0.85, 0.34])

    ax1 = timelabel(ax1, args.ttype)
    ax1.tick_params(axis="x", labelsize=12)
    ax1.set_ylabel("T mag [mag]")
    ax1.invert_yaxis()
    ax1.set_ylim([ymin1_l, ymax1_l])

    ax1_r = ax1.twinx()
    ax1_r.set_ylabel("Phase angle [deg]", color=col_a)
    ax1_r.set_ylim([ymin1_r, ymax1_r])
    # Line
    ax1_r.spines['right'].set_color(col_a)
    # Label numbers 
    ax1_r.tick_params(axis='y', colors=col_a, which='both')

    ax2 = timelabel(ax2, args.ttype)
    ax2.tick_params(axis="x", labelsize=12)
    ax2.set_ylabel("$\Delta$ [au]")
    ax2.set_xlabel(f"Year-Month")
    ax2.set_ylim([ymin2_l, ymax2_l])

    ax2_r = ax2.twinx()
    ax2_r.set_ylabel(f"$r_h$", color=col_a)
    ax2_r.set_ylim([ymin2_r, ymax2_r])
    # Line
    ax2_r.spines['right'].set_color(col_a)
    # Label numbers 
    ax2_r.tick_params(axis='y', colors=col_a, which='both')

    x, y = -0.05, 0.5
    ax1.yaxis.set_label_coords(x, y)

    # temporally create date list
    utclist = [
      datetime.datetime.strptime(i,'%Y-%b-%d %H:%M') for i in df["datetime_str"]]

    obj = df.at[0, "targetname"]
    label = "3I/ATLAS"

    # V-magnitude
    ax1.plot(
      utclist, df[key_mag],  color=col_V, ls=ls_V, lw=2, marker="", 
      label=label)
    # elevation
    ax1_r.plot(
      utclist, df["alpha"],  color=col_a, ls=ls_a, lw=2, marker="", 
      label=None)
    
    # delta
    ax2.plot(
      utclist, df["delta"],  color=col_V, ls=ls_V, lw=2, marker="", 
      label=label)
    # r_h
    ax2_r.plot(
      utclist, df["r"],  color=col_a, ls=ls_a, lw=2, marker="", 
      label=None)
    ax1.yaxis.set_major_locator(MultipleLocator(1.0))
    out = f"I_fig_ephem.{args.outtype}"
    out = os.path.join(args.outdir, out)
    plt.savefig(out, dpi=200)
