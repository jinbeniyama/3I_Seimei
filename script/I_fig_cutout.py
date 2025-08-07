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
from matplotlib.patches import FancyBboxPatch


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


def make_cutout_all(fi_list, center_list, band_list, out):
    """Make all cutout images.

    Parameters
    ----------
    fi_list : str
        list of fits file
    center_list : array-like
        list of x and y coordinates of the object
    band_list : array-like
        list of bands
    out : str
        output filename
    """

    sig_min, sig_max = 0.5, 5
    x_text, y_text = 0.13, 0.84

    # Full width in arcmin
    wi_arcmin = 1.00

    fig = plt.figure(figsize=(16, 20))

    nrows, ncols = 6, 6
    padding = 0.01 

    ax_width = (1.0 - padding * (ncols + 1)) / ncols
    ax_height = (1.0 - padding * (nrows + 1)) / nrows

    axes = []
    for i in range(nrows):
        for j in range(ncols):
            left = padding + j * (ax_width + padding)
            bottom = 1.0 - ((i + 1) * (ax_height + padding))
            ax = fig.add_axes([left, bottom, ax_width, ax_height])
            axes.append(ax)

    N_block = 12
    for n in range(N_block):
        axes_n = axes[3*n:3*(n+1)]
        
        fi3 = fi_list[n]
        center3 = center_list[n]
        band3 = band_list[n]

        for i, band in enumerate(band3):
            ax = axes_n[i]
            center = center3[i]
            fi = fi3[i]
            ax.axis("off")

    
            if band == "g":
                plt.rcParams["image.cmap"] = "Greens_r"
                col_text = "green"
            elif band == "r":
                plt.rcParams["image.cmap"] = "Reds_r"
                col_text = "red"
            elif band == "i":
                plt.rcParams["image.cmap"] = "PuRd_r"
                col_text = "magenta"
            elif band == "z":
                plt.rcParams["image.cmap"] = "Purples_r"
                col_text = "purple"

            ax.text(
                x_text, y_text, band, size=28, 
                color=col_text, transform=ax.transAxes,
                horizontalalignment="center",
                verticalalignment="center"
                )

            # Add box 
            box_width, box_height = 0.17, 0.17
            box_x = x_text - box_width / 2
            box_y = y_text - box_height / 2
            box = FancyBboxPatch(
                (box_x, box_y),
                box_width,
                box_height,
                boxstyle="round,pad=0.02",
                linewidth=2,
                facecolor="white",
                edgecolor="black",
                transform=ax.transAxes,
            )
            ax.add_patch(box)

            hdu = fits.open(fi)
            hdr = hdu[0].header
            utc = hdr["UTC"].replace("T", " ")[:19]
            data = hdu[0].data
            ny, nx = data.shape

            # Set title in on the center panel
            if i == 1:
                ax.set_title(utc)

            # Background ======================================================
            data = data.byteswap().newbyteorder()

            bg_engine = sep.Background(data)
            bg_global = bg_engine.globalback
            # Used to estimate flux error
            bg_rms = bg_engine.globalrms

            # Just remove typical background
            data -= bg_global
            # Background ======================================================

            xc, yc = center
            print(f"    Original center (xc, yc) = ({xc:.1f}, {yc:.1f})")
            # Centroid
            sigma = 2
            data = data.astype(data.dtype.newbyteorder('='))
            xc, yc, wflag = sep.winpos(data, xc, yc, sigma)
            print(f"    Winpos center (xc, yc) = ({xc:.1f}, {yc:.1f})")


            # Cut the data ====================================================
            # 2. by xcenter, ycenter, and width
            # When (xc, yc) = (51, 51) and width = 50,
            # cut region is in 1 <= x <= 101 and 1 <= y <= 101 in pixel.
            # Cut region is in 0 <= x <= 100 and 0 <= y <= 100 in python.
            # width should be even

            # Harf width in pixel (0.35 arcsec/pix)
            wi = int(wi_arcmin*60/0.35/2.)

            xmin, xmax = xc - wi - 1, xc + wi
            ymin, ymax = yc - wi - 1, yc + wi
            xmin, xmax = int(xmin), int(xmax)
            ymin, ymax = int(ymin), int(ymax)
            print(f"Data Region (x, y)=({xmin}-{xmax}, {ymin}-{ymax})")
            print(f"Data length (Nx, Ny)=({xmax-xmin}, {ymax-ymin})")

            data_cut = data[ymin:ymax, xmin:xmax]
            # Cut the data ====================================================

            # Remove local background again
            bg_global_local = np.median(data_cut)
            data -= bg_global_local

            # Flip ud (North is up and East is left in a figure.)
            #data = np.flipud(data)
            err = bg_rms
            print(
              f"(mean, std)=({np.mean(data):.1f}, {np.std(data):.1f})"
              )

            # Plot ============================================================
            vmin, vmax = sig_min*bg_rms, sig_max*bg_rms
            im = ax.imshow(data, vmin=vmin, vmax=vmax)

            # Circle
            rad_circle = 10
            circle = Circle(
                (xc, yc), rad_circle,
                edgecolor='white', facecolor='none', linewidth=1.5)
            ax.add_patch(circle)

            # Check range
            ax.set_xlim([xmin, xmax])
            ax.set_ylim([ymin, ymax])


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


    # Show all cutout images
    #  1. gri1
    #  2. grz1
    #  3. gri2
    #  4. grz2
    #  5. gri3
    #  6. grz3
    #  7. grz4
    #  8. grz5
    #  9. gri4
    # 10. gri5
    # 11. gri6
    # 12. gri7
    
    # Save fi_list (fits), co_list (coordinate of 3I), and band_list (bands)
    stack = "mean"
    fi_list = []
    co_list = []
    band_list = []
    # 2025-07-15 g, r, i, 600 s, 7 set
    co_gri_list = [
        [(1080, 643), (1066, 653), (1078, 633)],
        [(1082, 635), (1068, 645), (1079, 625)],
        [(1081, 635), (1067, 645), (1077, 625)],
        [(1078, 630), (1064, 640), (1075, 620)],
        [(1078, 637), (1064, 647), (1075, 627)],
        [(1082, 640), (1068, 648), (1079, 629)],
        [(1085, 630), (1070, 640), (1085, 620)],
    ]
    fi_lists = []
    for i in range(len(co_gri_list)):
        f_g = f"/Users/beniyama/research/Seimei_3I/Seimei202507_3I/20250715/3I_gri/shift_600s_{stack}/nonsid_3I_gri_600s_0_00{i+1}.fits"
        f_r = f"/Users/beniyama/research/Seimei_3I/Seimei202507_3I/20250715/3I_gri/shift_600s_{stack}/nonsid_3I_gri_600s_1_00{i+1}.fits"
        f_i = f"/Users/beniyama/research/Seimei_3I/Seimei202507_3I/20250715/3I_gri/shift_600s_{stack}/nonsid_3I_gri_600s_2_00{i+1}.fits"
        fi_list.append([f_g, f_r, f_i])
        band_list.append(["g", "r", "i"])
        co_list.append(co_gri_list[i])
    
    # 2025-07-15 g, r, z, 600 s, 5 set
    co_grz_list = [
        [(1080, 643), (1066, 653), (1078, 633)],
        [(1088, 643), (1074, 652), (1086, 632)],
        [(1085, 642), (1070, 652), (1082, 635)],
        [(1076, 633), (1062, 643), (1073, 623)],
        [(1076, 638), (1062, 647), (1073, 627)],
        ]
    
    for i in range(len(co_grz_list)):
        f_g = f"/Users/beniyama/research/Seimei_3I/Seimei202507_3I/20250715/3I_grz/shift_600s_{stack}/nonsid_3I_grz_600s_0_00{i+1}.fits"
        f_r = f"/Users/beniyama/research/Seimei_3I/Seimei202507_3I/20250715/3I_grz/shift_600s_{stack}/nonsid_3I_grz_600s_1_00{i+1}.fits"
        f_z = f"/Users/beniyama/research/Seimei_3I/Seimei202507_3I/20250715/3I_grz/shift_600s_{stack}/nonsid_3I_grz_600s_2_00{i+1}.fits"
        fi_list.append([f_g, f_r, f_z])
        band_list.append(["g", "r", "z"])
        co_list.append(co_grz_list[i])

    
    N_fi =  len(fi_list)
    N_co =  len(co_list)
    N_band =  len(band_list)

    # Sort by obs. time
    order_obs = [
        1, 3, 5, 9, 10, 11, 12, 
        2, 4, 6, 7, 8]

    fi_list_sort = [None] * len(fi_list)
    co_list_sort = [None] * len(co_list)
    band_list_sort = [None] * len(band_list)
    for i, pos in enumerate(order_obs):
        fi_list_sort[pos - 1] = fi_list[i]
        co_list_sort[pos - 1] = co_list[i]
        band_list_sort[pos - 1] = band_list[i]

    # Plot all 
    # 6 x 6
    out = f"I_fig_cutout_all.{args.outtype}"
    out = os.path.join(outdir, out)
    make_cutout_all(
        fi_list_sort, co_list_sort, band_list_sort, out)
