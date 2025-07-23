#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
!!! This script is for JB's MacBook. !!!

Create cutout images of 3I/ATLAS on 2025-07-15 gri/z images.
"""
import os 
from argparse import ArgumentParser as ap
import numpy as np
import matplotlib.pyplot as plt  
from astropy.io import fits as fits
import sep
from matplotlib.patches import Circle


def make_cutout(fi, center, band):
    """Make cutout image.

    Parameters
    ----------
    fi : str
        fits file
    center : array-like
        x and y coordinates of the object
    band : str
        "g", "r", "i", or "z"
    """
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_axes([0.00, 0.00, 1.0, 1.0])

    sig_min, sig_max = 0.5, 5
    ax.axis("off")
    
    x_text, y_text = 0.15, 0.85
    if band == "g":
        plt.rcParams["image.cmap"] = "Greens_r"
    elif band == "r":
        plt.rcParams["image.cmap"] = "Reds_r"
    elif band == "i":
        plt.rcParams["image.cmap"] = "PuRd_r"
    elif band == "z":
        plt.rcParams["image.cmap"] = "Purples_r"
    ax.text(
        x_text, y_text, band, size=30, 
        color="white", transform=ax.transAxes)

    hdu = fits.open(fi)
    hdr = hdu[0].header
    data = hdu[0].data
    assert len(data.shape)==2, "Input 2-d fits!"
    ny, nx = data.shape

    # Cut the data ========================================================
    # 2. by xcenter, ycenter, and width
    # When (xc, yc) = (51, 51) and width = 50,
    # cut region is in 1 <= x <= 101 and 1 <= y <= 101 in pixel.
    # Cut region is in 0 <= x <= 100 and 0 <= y <= 100 in python.
    # width should be even

    # Full width in arcmin
    wi_arcmin = 1.00
    # Harf width in pixel (0.35 arcsec/pix)
    wi = wi_arcmin*60./0.35/2.
    xc, yc = center
    xmin, xmax = xc - wi - 1, xc + wi
    ymin, ymax = yc - wi - 1, yc + wi
    xmin, xmax = int(xmin), int(xmax)
    ymin, ymax = int(ymin), int(ymax)
    print(f"Data Region (x, y)=({xmin}-{xmax}, {ymin}-{ymax})")
    print(f"Data length (Nx, Ny)=({xmax-xmin}, {ymax-ymin})")

    data = data[ymin:ymax, xmin:xmax]
    # Calculate aspect ratio
    aspect_ratio = (xmax-xmin)/(ymax-ymin)
    print(f"aspect ratio is {aspect_ratio:.4f}")
    # Cut the data ========================================================


    # Remove background ===================================================
    data = data.byteswap().newbyteorder()
    bg_engine = sep.Background(data)
    bg_engine.subfrom(data)
    bg_global = bg_engine.globalback
    bg_rms = bg_engine.globalrms
    bg_info = {'level': bg_global, 'rms': bg_rms}
    data -= np.median(data)

    # Flip ud (North is up and East is left in a figure.)
    data = np.flipud(data)
    err = np.round(bg_info["rms"],2)
    print(
      f"(mean, std)=({np.mean(data):.1f}, {np.std(data):.1f})"
      )
    # Remove background ===================================================

      
    # Plot ================================================================
    vmin, vmax = sig_min*err, sig_max*err
    im = ax.imshow(data, vmin=vmin, vmax=vmax)

    # Circle
    rad_circle = 6
    circle = Circle(
        (wi, wi), rad_circle,
        edgecolor='white', facecolor='none', linewidth=2)
    ax.add_patch(circle)

    out = f"I_fig_cutout_{band}.{args.outtype}"
    out = os.path.join(outdir, out)
    plt.savefig(out)


if __name__ == "__main__":
    parser = ap(
        description="Create cutout of 3I/ATLAS")
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
    
    # Fits files (JB's MBP)
    # g,r,i,z
    bands = ["g", "r", "i", "z"]
    # 60 s x 10 = 600 s
    flist = [
        "/Users/beniyama/research/Seimei_3I/Seimei202507_3I/20250715/3I_gri/shift_600s_mean/nonsid_3I_gri_600s_0_004.fits",
        "/Users/beniyama/research/Seimei_3I/Seimei202507_3I/20250715/3I_gri/shift_600s_mean/nonsid_3I_gri_600s_1_004.fits",
        "/Users/beniyama/research/Seimei_3I/Seimei202507_3I/20250715/3I_gri/shift_600s_mean/nonsid_3I_gri_600s_2_004.fits",
        "/Users/beniyama/research/Seimei_3I/Seimei202507_3I/20250715/3I_grz/shift_600s_mean/nonsid_3I_grz_600s_2_002.fits"]
    # central positions with JB's own eyes
    center = [
        (1078, 630), (1064, 640), (1075, 620), (1086, 632)
        ]

    for (fi, cen, band)  in zip(flist, center, bands):
        print(f"{fi}, {cen}, {band}")
        make_cutout(fi, cen, band)
