#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Useful info. and functions for 3I paper.
"""
import os
import pandas as pd
import numpy as np
import datetime
from decimal import Decimal, ROUND_HALF_UP
from scipy.optimize import curve_fit
from astroquery.jplhorizons import Horizons
import matplotlib.pyplot as plt  
from astropy import units as u
from astropy.time import Time
from matplotlib.dates import DateFormatter, HourLocator, DayLocator, MonthLocator
import classy
from matplotlib.patches import Rectangle


# Constants and our results ===================================================
loc_Seimei = {
    "lon":        133.5967, 
    "lat":         34.5769, 
    "elevation" :   0.355
    }

# Meech+2017, SDSS (not used)
I1_colors_Meech = dict(
      g_r = 0.84, g_rerr = 0.05,
      # g-i = 1.15+-0.10
      # g-z = 1.25+-0.10
      # -> r-i = g-i - (g-r) = 0.31 (+-0.11)
      # -> i-z = g-z - (g-i) = 0.10 (+-0.14)
      # -> r-z = g-z - (g-r) = 0.41 (+-0.11)
      r_i = 0.31, r_ierr = 0.11,
      i_z = 0.10, i_zerr = 0.14,
      r_z = 0.41, r_zerr = 0.11,
    )

# Seligman+2025, Pan-STARRS
# July 2 and 4, 2025, 
# Faulkes Telescope North(FTN)/MuSCAT3 
# (FTS was used for lightcurve obs.)
# aperture radii 1.46 arcsec
I_colors_Seligman = dict(
      g_r = 0.85, g_rerr = 0.03,
      r_i = 0.25, r_ierr = 0.03,
      i_z = 0.20, i_zerr = 0.08,
      # From r–i=0.25±0.03, i–z=–0.20±0.08
      r_z = 0.45, r_zerr = 0.085,
    )

# Bolin+2025, SDSS, 
#   KAO 1.88 m      UBVR 2025-07-02
#   200 inch/NGPS   gri  2025-07-03
#   3.5 m/ARCTIC    griz 2025-07-06
# Note: aperture radii 4 arcsec (10,000 km from
#       the comet’s nucleus at its 3.36 - 3.45 au distance from the Earth
I_colors_Bolin = dict(
      g_r = 0.84, g_rerr = 0.05,
      r_i = 0.16, r_ierr = 0.03,
      i_z = -0.02, i_zerr = 0.07,
      # From r–i=0.16±0.03, i–z=–0.02±0.07
      r_z = 0.14, r_zerr = 0.076,
    )
# 3.5 m/ARCTIC
I_mags_Bolin = dict(
    g = 18.58, gerr = 0.06,
    r = 17.77, rerr = 0.02,
    i = 17.59, ierr = 0.04,
    z = 17.61, zerr = 0.05,
    )
# 200 inch/NGPS
# g=18.72+-0.02
# r=17.86+-0.01
# i=17.72+-0.02

# Kareta+2025, SDSS
#  IRTF r 2025-07-03
#  IRTF gi 2025-07-04
# Note: 4 arcsec radius circular aperture
# They don't believe the mean values of I_mags_Kareta are correct
# since the field was crowded.
# So only g-i is reliable I think.
I_mags_Kareta = dict(
    g = 18.39, gerr = 0.02,
    r = 17.74, rerr = 0.01,
    i = 17.37, ierr = 0.02,
    )
# From best frames 
# PS (written in the summary)
#   g-i = 1.16 +- 0.20 (main text and summary)
#   (typo? g-i = 1.06 +- 0.11 (abstract, main text))
I_colors_Kareta = dict(
    g_i = 1.16, g_ierr = 0.20,
    )

# Puzia+2025, PS (SDSS 0.86+-0.05, PS=0.73+-0.05)
# 4.1 m SOAR Telescop/Goodman High Throughput Spectrograph
# 2025-07-04 UT 04:50:33.0 to 06:36:01.0
# JB confirmed that these colors are consistent with SDSS2PS_col function.
# Estimated from spectrum
I_colors_Puzia = dict(
      # SDSS
      #g_r = 0.86, g_rerr = 0.05,
      # PS
      g_r = 0.73, g_rerr = 0.05,
    )


# Santana-Ros+2025, 
# griz colors are calibrated to ATLAS All-Sky Stellar Reference Catalog 
# (Tonry et al. 2018, AJ, private communication w/Toni),
# in which all griz photometry has been transformed to the 
# Pan-STARRS gP1, rP1, iP1, and zP1 bandpasses.
# Refcat2 magnitudes are effectively Pan-STARRS DR1 magnitudes 
# (north of decl. -30 deg and mag > 14).
# (see 7. Summary and Conclusions of Tonry et al., 2018, AJ)
#   1. Faulkes Telescope North (FTN)/MuSCAT3, griz
#   2. Faulkes Telescope South (FTS)/MuSCAT4, griz
#   3. 1-meter Lesedi telescope, ugriz
#   (4. Skalnaté Pleso Observatory, V-R)
# 2025-07-02 to 2025-07-29
# Avaraged values
# Note: aperture radii are various.
I_colors_Toni = dict(
      g_r = 0.65, g_rerr = 0.03,
      r_i = 0.27, r_ierr = 0.03,
      i_z = 0.10, i_zerr = 0.04,
      # From r–i=0.27±0.03, i–z=0.10±0.04
      # -> Calculate after conversion
      #r_z = 0.37, r_zerr = 0.05,
    )

def lam_W18(band):
    """Return wavelength and its fwhm of a filter.

    Parameter
    ---------
    band : str
        filter

    Returns
    -------
    lam_mean, lam_fwhm
        pivot wavelfnth and its fwhm
    """
    # Willmer 2018, ApJS, 236, 47. Table 4
    # lambda pivot, fwhm
    # e.g, Extract mean, minimum, and maximum wavelength
    # $ lam, lam0, lam1, = lam_W18["g"]
    # Johnson V
    V_lam, V_fwhm = 0.5511, 0.0812
    ## Pan-SRTARRS 1
    g_lam, g_fwhm = 0.4849, 0.1256
    r_lam, r_fwhm = 0.6201, 0.1404
    i_lam, i_fwhm = 0.7535, 0.0698
    z_lam, z_fwhm = 0.8674, 0.1034
    Y_lam, Y_fwhm = 0.9628, 0.0629
    ## 2MASS
    J_lam,  J_fwhm  = 1.2393, 0.2027
    H_lam,  H_fwhm  = 1.6495, 0.2610
    Ks_lam, Ks_fwhm = 2.1638, 0.2785

    if band == "V":
        lam_mean, lam_fwhm = V_lam, V_fwhm
    elif band == "g":
        lam_mean, lam_fwhm = g_lam, g_fwhm
    elif band == "r":
        lam_mean, lam_fwhm = r_lam, r_fwhm
    elif band == "i":
        lam_mean, lam_fwhm = i_lam, i_fwhm
    elif band == "z":
        lam_mean, lam_fwhm = z_lam, z_fwhm
    elif band == "Y_PS":
        lam_mean, lam_fwhm = Y_lam, Y_fwhm
    elif band == "J":
        lam_mean, lam_fwhm = J_lam, J_fwhm
    elif band == "H":
        lam_mean, lam_fwhm = H_lam, H_fwhm
    elif band == "Ks":
        lam_mean, lam_fwhm = Ks_lam, Ks_fwhm
    else:
        assert False, "Not implemented."

    return lam_mean, lam_fwhm


def sun_mag_W18(band, system):
    """Return absolute magnitude of the Sun.

    Parameters
    ----------
    band : str
        filter
    system : str
        vega or AB

    Return
    ------
    sun_mag : float 
        absolute magnitude of the Sun
    """
    # Solar color (Table 3 of Willmer 2018)
    # V  : Johnson V
    # Y  : PS1
    # JHK: 2MASS
    
    if system == "vega":
        sun_mags = dict(
            u_SDSS=5.49,
            g_SDSS=5.23,
            r_SDSS=4.53,
            i_SDSS=4.19,
            z_SDSS=4.01,
            g_PS=5.14,
            r_PS=4.53,
            i_PS=4.18,
            z_PS=4.02,
            Y_PS=3.99,
            V=4.81,
            J=3.67,
            H=3.32,
            Ks=3.27,
            )
    elif system == "AB":
        sun_mags = dict(
            u_SDSS=6.39,
            g_SDSS=5.11,
            r_SDSS=4.65,
            i_SDSS=4.53,
            z_SDSS=4.50,
            g_PS=5.03,
            r_PS=4.64,
            i_PS=4.52,
            z_PS=4.51,
            Y_PS=4.50,
            V=4.80,
            J=4.54,
            H=4.66,
            Ks=5.08,
            )

    sun_mag = sun_mags[band]
    return sun_mag
# Constants and our results ===================================================


def add_obsinfo(df, loc):
    """Add obs. info.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with 
         'obj', 'obsdate0', 'obsdate1', 'code', 
    loc : str or dict
        location for JPL query

    Return
    ------
    df : pandas.DataFrame
        'H', 'V0', 'V1', 'elev0', 'elev1', 'vnorm0', 'vnorm1', 'r0', 'r1',
        'delta0', 'delta1', 'alpha0', 'alpha1' 
        added 15-columns DataFrame
    """

    columns = df.columns.tolist()
    assert "obj" in columns, "Check input."
    assert "obsdate0" in columns, "Check input."

    
    df = df.assign(
        H=0.0, V0=0.0, V1=0.0, elev0=0.0, vnorm0=0.0,
        r0=0.0, r1=0.0, delta0=0.0, alpha0=0.0, lunarillum0=0.0, lunarelong0=0.0,
        glon=0.0, glat=0.0)
    for idx, row in df.iterrows():
        print(f"  Query {idx+1:02d} start")
        obj = row["obj"]
        jd0 = row["jd0"]
        obj = f"{obj}"

        # Query with a single time
        jpl = Horizons(id=obj, location=loc, epochs=jd0)

        eph= jpl.ephemerides()
        # Visible magnitude
        df.at[idx, "V0"] = eph["Tmag"][0]
        # Elevation
        df.at[idx, "elev0"] = eph["EL"][0]
        df.at[idx, "elev1"] = eph["EL"][-1]
        # minimum/maximum airmass (do not corresponds to first and last frames)
        airmass = 1/np.cos(np.radians(90-eph["EL"]))
        df.at[idx, "airmass0"] = np.min(airmass)
        dcos = np.cos(np.radians(eph["DEC"][0]))

        # Already arcsec/hour in sky motion. Do not need cos correction
        vnorm0 = (eph["RA_rate"][0]**2 + eph["DEC_rate"][0]**2)**0.5
        # arcsec/hour to arcsec/s
        df.at[idx, "vnorm0"] = vnorm0/3600.

        #  'PsAng   PsAMV' =
        #      The position angles of the extended Sun-to-target radius vector ("PsAng")
        #  and the negative of the targets' heliocentric velocity vector ("PsAMV"), as
        #  seen in the observers' plane-of-sky, measured counter-clockwise (east) from
        #  reference-frame north-pole. Primarily intended for ACTIVE COMETS, "PsAng"
        #  is an indicator of the comets' gas-tail orientation in the sky (being in the
        #  anti-sunward direction) while "PsAMV" is an indicator of dust-tail orientation.
        #  Units: DEGREES
        df.at[idx, "phi0"] = eph["sunTargetPA"][0]

        # Heliocentric distance
        df.at[idx, "r0"] = eph["r"][0]
        # Geocentric distance
        df.at[idx, "delta0"] = eph["delta"][0]
        # Solar phase angle
        df.at[idx, "alpha0"] = eph["alpha"][0]

        # Percent of the Moon Illuminated and elongation 
        df.at[idx, "lunarelong0"] = eph["lunar_elong"][0]
        df.at[idx, "lunarillum0"] = eph["lunar_illum"][0]

        # Galactic longitude and latitude
        df.at[idx, "glon0"] = eph["GlxLon"][0]
        df.at[idx, "glat0"] = eph["GlxLat"][0]

        
    return df


def tab_phot(df, date, out):
    """Create obs information tex table from 15-columns DataFrame.

    Parameters
    ----------
    df : pandas.DataFrame
      15 columns NEOs DataFrame
    date : str
        date like 2024-06-16
    out : str
      output filename
    """
    
    # Like 2023 Mar 24
    # Time like 2022 Dec 20 15:30:37--16:29:49
    date = datetime.datetime.strptime(date, "%Y-%m-%d")
    # Time like 2022 December 05
    date = datetime.datetime.strftime(date, "%B %d, %Y")

    note = (
            r"""
    Observation time in UT in midtime of exposure (Obs. Date), 
    and filters (Filter) are listed.
    Predicted V band apparent magnitudes (V)
    at the observation starting time
    are referred to NASA Jet Propulsion Laboratory (JPL) Horizons
    """
    +
    f" as of {date}."
    +
    r"""
    Elevations of 3I/ATLAS to calculate air mass range (Air Mass) are 
    also referred to NASA JPL Horizons.
            """)

    head = (
        r"""\begin{longtable}{rcccl}
\caption{Summary of the observations}\label{tab:obs}
\hline\noalign{\vskip3pt}
        Obs. Date & Filter & V     & Air Mass   & Note \\ [2pt]
        (UTC)     &        & (mag) &            & \\ [2pt]
\hline\noalign{\vskip3pt}
\endfirsthead
\hline\noalign{\vskip3pt}
        Obs. Date & Filter & V     & Air Mass   & Note \\ [2pt]
        (UTC)     &        & (mag) &            & \\ [2pt]
\hline\noalign{\vskip3pt}
\endhead
\hline\noalign{\vskip3pt}
\endfoot
\hline\noalign{\vskip3pt}
\multicolumn{2}{@{}l@{}}{\hbox to0pt{\parbox{160mm}{\footnotesize
\hangindent6pt\noindent
\hbox to6pt{\footnotemark[$*$]\hss}\unskip%"""
+ note +
"""}\hss}}\endlastfoot
""")

    with open(out, "w") as f:
        f.write(head)
        for idx, row in df.iterrows():
            obj_tex = row["obj"]
            
            # Not included in the table
            r = row["r0"]
            delta = row["delta0"]
            alpha = row["alpha0"]
            lunarillum = row["lunarillum0"]
            lunarelong = row["lunarelong0"]
            glon       = row["glon0"]
            glat       = row["glat0"]
            # arcsec/s
            vnorm      = row["vnorm0"]
            print(f"  {idx+1}: Heliocentric distance  {r:.2f}")
            print(f"           Observer-centric dist. {delta:.2f}")
            print(f"           Phase angle            {alpha:.2f}")
            print(f"           Lunar illumination     {lunarillum:.2f} %")
            print(f"           Lunar elongation       {lunarelong:.2f} deg")
            print(f"           Galactic lon, lat      ({glon:.2f}, {glat:.2f})")
            print(f"           Sky motion             {vnorm:.2f} arcsec/s")
            print(f"                                  {vnorm*60:.2f} arcsec/min")
            
            # Filter 
            fltr = row["band"]
            # gri or grz
            if len(fltr)==3:
                b1, b2, b3 = fltr[0], fltr[1], fltr[2]
                fltr_str = f"${b1},{b2},{b3}$"
            else:
                fltr_str = fltr

            # Time like 2022 Dec 20 15:30:37--16:29:49
            t0 = row["obsdate0"].replace("T", " ")[:19]
            # 2022-12-20
            t0_head = t0[0:10]
            t0_head = datetime.datetime.strptime(t0_head, "%Y-%m-%d")
            # 2022 Dec 20
            t0_head = datetime.datetime.strftime(t0_head, "%Y %b %d")
            t0_tail = t0[11:19]

            obstimeinfo = f"{t0_head} {t0_tail}"
            # Do not show twice
            if idx > 0:
                if t0_head == t0_head0:
                    obstimeinfo = f"{t0_tail}"
                else:
                    pass

            text = (
                f"{obstimeinfo}&"
                f" {fltr_str} & {row['V0']:.1f} &"
                f" {row['airmass0']:.2f} &"
                f" {row['note']}"
                r"\\")
            f.write(text)
            f.write("\n")
            t0_head0 = t0_head
        f.write(r"\end{longtable}")
# To make table ===============================================================


# Plot ========================================================================
mycolor = [
    "#AD002D", "#1e50a2", "#69821b", "#f055f0", "#afafb0", 
    "#0095b9", "#89c3eb", "#ec6800", "cyan", "gold",
    "magenta"
    ] 
mycolor = mycolor*500

# linestyle
myls = ["solid", "dashed", "dashdot", "dotted", 
        (0, (5, 3, 1, 3, 1, 3)), (0, (4,2,1,2,1,2,1,2)),
        (0, (4,2,1,2,1,2,1,2,1,2)), (0, (4,2,1,2,1,2,1,2,1,2,1,2)),
        ]
myls = myls*100

# marker
mymark = ["o", "^", "s", "D", "*", "v", "<", ">", "h", "x"]
mymark = mymark*500

bandcolor = {
    "g": "#69821b", "r": "#AD002D", "i": "magenta", "z": "#9400d3"}
bandmark = {
    "g": mymark[0], "r": mymark[1], "i": mymark[2], "z": mymark[3]}


def mygrid(ax):
    ax.grid(ls="dashed", lw=0.4)

def mymathtext():
    """Update config of matplotlib for beautiful text.
    """
    from matplotlib import _mathtext as mathtext
    mathtext.FontConstantsBase = mathtext.ComputerModernFontConstants
    # Space after the sub/superscript
    mathtext.FontConstantsBase.script_space = 0.01
    # Space between the str and sub/superscript
    mathtext.FontConstantsBase.delta = 0.01
    mathtext.FontConstantsBase.subdrop = 0
    # The same baseline for subscript
    mathtext.FontConstantsBase.sub1 = 0
    # Scale of sub/superscripts
    mathtext.SHRINK_FACTOR = 0.6


def stype2colmark(stype):
    """Return marker for input spectral type.

    Parameter
    ---------
    stype : str
        spectral type

    Returns
    -------
    color : str
        color
    marker : str
        marker
    """
    if stype=="S":
        color = mycolor[0]
        mark = mymark[0]
    elif stype=="V":
        color = mycolor[3]
        mark = mymark[1]
    elif stype=="X":
        color = mycolor[2]
        mark = mymark[2]
    elif stype=="K":
        color = mycolor[7]
        mark = mymark[3]
    elif stype=="L":
        color = mycolor[9]
        mark = mymark[5]
    elif stype=="C":
        color = mycolor[1]
        mark = mymark[6]
    elif stype=="B":
        color = mycolor[5]
        mark = mymark[7]
    elif stype=="D":
        color = mycolor[6]
        mark = "p"
    elif stype=="A":
        color = mycolor[4]
        mark = mymark[8]
    # No Q-type in SDSSMOC
    else:
        color, mark = "black", "o"
    return color, mark


# Plot ========================================================================


# Error =======================================================================
def log10err(val, err):
    """Calculate log10 error.
    """
    return err/val/np.log(10)


def diverr_series(val1, err1, val2, err2):
    """Calculate error for division.

    Parameters
    ----------
    val1 : float
        value 1
    err1 : float
        error 1
    val2 : float
        value 2
    err2 : float
        error 2

    Return
    ------
    err : float
        propageted error
    """
    err_list = []
    for v1, e1, v2, e2 in zip(val1, err1, val2, err2):
        err = np.sqrt((e1*1.0/v2)**2 + (e2*v1/v2**2)**2)
        err_list.append(err)
    return err_list


def mulerr(val1, err1, val2, err2):
    """Calculate error for multiple.

    Parameters
    ----------
    val1 : float
        value 1
    err1 : float
        error 1
    val2 : float
        value 2
    err2 : float
        error 2
    """
    return np.sqrt((val2*err1)**2 + (val1*err2)**2)


def adderr(*args):
    """Calculate additional error.

    Parameters
    ----------
    args : array-like
        list of values

    Return
    ------
    err : float
        calculated error
    """
    err = np.sqrt(np.sum(np.square(args)))
    return err


def adderr_series(*args):
    """Add error of multiple pandas.Series.

    Parameters
    ----------
    args : array-like
        list of pandas.Series 

    Return
    ------
    err_s : pandas.Series
        single pandas.Series of calculated errors
    """ 
    for i,x in enumerate(args):
        assert type(x)==type(pd.Series()), "Sould be Series"
        if i==0:
            temp = x.map(np.square)
        else:
            temp += x.map(np.square)
    err_s = temp.map(np.sqrt)
    return err_s
# Error =======================================================================


# Colors ======================================================================
# Tonry2012
def gr2V_mag(g, gerr, r, rerr):
    """Convert g and r mag to V-band magnitude based on Tonry+2012.

    Parameters
    ----------
    g : float
        g-band magnitude
    gerr : float
        g-band magnitude error
    r : float
        r-band magnitude
    rerr : float
        r-band magnitude error

    Returns
    -------
    V : float
        V-band magnitude
    Verr : float
        V-band magnitude error
    """
    # Equation (2) in block 4 in table 6
    V = r + 0.006 + 0.474*(g-r)
    # Error of colors
    grerr = adderr(gerr, rerr)
    # rerr, error of conversion, and error of colors
    Verr = adderr(rerr, 0.012, 0.474*grerr)
    return V, Verr


def SDSS2PS_mag(df, key0=None, key1=None):
    """Convert SDSS magnitude to Pan-STARRS magnitude.

    (convert origianl 'g' to 'g_SDSS', and add converted magnitude as 'g_PS')
    See block 1 in table 6 in Tonry+2012.
    Note: 'err's are error for conversion itself.
           ex) when g1 = g0 + (g0-r0)*C, and conversion error is err, 
               total error (err_total) is 
                 err_total**2 = err**2 + (C*g0err)**2 + (C*r0err)**2

    Parameters
    ----------
    df : pandas.DataFrame
        dataframe contains SDSS g,r,i,z-band magnitude
    key0 : dict, optional
        magnitude and magnitude error keys of original DataFrame
    key1 : dict, optional
        magnitude and magnitude error keys to be added

    Return
    ------
    df : pandas.DataFrame
        dataframe contains Pan-STARRS g,r,i,z-band magnitude
    """
    # Data dimension
    n = len(df)

    # Original keys to be converted.
    if key0 is None:
      key0 = {
        "g":"g", "gerr":"gerr",
        "r":"r", "rerr":"rerr",
        "i":"i", "ierr":"ierr",
        "z":"z", "zerr":"zerr"
        }
    # New keys.
    if key1 is None:
      key1 = {
        "g":"g_PS", "gerr":"gerr_PS",
        "r":"r_PS", "rerr":"rerr_PS",
        "i":"i_PS", "ierr":"ierr_PS",
        "z":"z_PS", "zerr":"zerr_PS"
       }

    # Estimate g, r, i, and z
    # Coefficients from Table 6 in Tonry+2012
    coeff = {
        "B0":  [-0.012,  0.000,  0.004, -0.013], 
        "B1":  [-0.139, -0.007, -0.014,  0.039], 
        "Berr": [ 0.007,  0.002,  0.003,  0.009]
        }
    bands = ["g", "r", "i", "z"]
    for idx, b in enumerate(bands):
        B0 = coeff["B0"][idx]
        B1 = coeff["B1"][idx]
        err = coeff["Berr"][idx]
        # Use g-r color for all conversion (see Table 6)
        df[key1[b]] = df[key0[b]] + B0 + (df[key0["g"]] - df[key0["r"]])*B1
        print(f"  Conversion: {key1[b]} = {key0[b]} + {B0} + (g-r)*{B1}")
        ## Create a pandas.Series of a conversion error
        converr = [err]*n
        converr_series = pd.Series(data=converr)
        df[key1[f"{b}err"]] = adderr_series(
            df[key0[f"{b}err"]], converr_series, 
            B1*df[key0["gerr"]], B1*df[key0["rerr"]])
    return df


def SDSS2PS_col(df, key0=None, key1=None):
    """
    Convert colors in the SDSS to Pan-STARRS.
    See block 1 in table 6 in Tonry+2012.
    Note: 'err's are error for conversion itself.
          When 
          Color1 = C0 + C1*(Color0)*C , and conversion error is err1 and err2
            (two errors since the equations are some of two eqs. in Table 6)
          err_total**2 = err1**2 + err2**2 + (C*(g0-r0)err)**2

    Parameter
    ---------
    df : pandas.DataFrame
        dataframe contains gri color
    key0 : dict, optional
        color and color error keys of original DataFrame
    key1 : dict, optional
        color and color error keys to be added

    Return
    ------
    df : pandas.DataFrame
        dataframe contains PS color
    """
    # Data dimension
    n = len(df)

    # Original keys to be converted.
    if key0 is None:
        key0 = {
            "g_r":"g_r_SDSS", "g_rerr":"g_rerr_SDSS",
            "r_i":"r_i_SDSS", "r_ierr":"r_ierr_SDSS",
            "i_z":"i_z_SDSS", "i_zerr":"i_zerr_SDSS",
        }
    # New keys.
    if key1 is None:
        key1 = {
            "g_r":"g_r_PS", "g_rerr":"g_rerr_PS",
            "r_i":"r_i_PS", "r_ierr":"r_ierr_PS",
            "i_z":"i_z_PS", "i_zerr":"i_zerr_PS",
        }

    # Coefficients from block 1 in Table 6 (Tonry+2012)
    coeff = {
        "B0":   [ -0.012,  0.000,   0.004, -0.013], 
        "B1":   [ -0.139, -0.007,  -0.014,  0.039], 
        "Berr": [  0.007,  0.002,   0.003,  0.009]
        }

    # Calculate (g-r)PS from equations (1) - (2)
    df[key1["g_r"]] = (
        (coeff["B0"][0] - coeff["B0"][1]) 
        + (coeff["B1"][0] - coeff["B1"][1] + 1)*df[key0["g_r"]]
        )
    # Sum of err[0], err[1], and (B1[0]-B1[1]+1)*g_rerr
    err1 = [coeff["Berr"][0]]*n
    err2 = [coeff["Berr"][1]]*n
    err1_s = pd.Series(data=err1)
    err2_s = pd.Series(data=err2)
    err_gr = df[key0["g_rerr"]]
    df[key1["g_rerr"]] = adderr_series(
        err1_s, err2_s, (coeff["B1"][0]-coeff["B1"][1]+1)*err_gr)


    # Calculate (r_i)PS from equations (2) - (3)
    df[key1["r_i"]] = (
      (coeff["B0"][1] - coeff["B0"][2]) 
      + (coeff["B1"][1] - coeff["B1"][2])*df[key0["g_r"]] + df[key0["r_i"]]
      )
    # Sum of err[1], err[2], (B1[1]-B1[2])*g_rerr, and r_ierr
    err3 = [coeff["Berr"][2]]*n
    err3_s = pd.Series(data=err3)
    err_ri = df[key0["r_ierr"]]
    df[key1["r_ierr"]] = adderr_series(
      err2_s, err3_s, (coeff["B1"][1]-coeff["B1"][2])*err_gr, err_ri)


    # Calculate (i_z)PS from equations (3) - (4)
    df[key1["i_z"]] = (
      (coeff["B0"][1] - coeff["B0"][2]) 
      + (coeff["B1"][1] - coeff["B1"][2])*df[key0["g_r"]] + df[key0["i_z"]]
      )
    # Sum of err[1], err[2], (B1[1]-B1[2])*g_rerr, and i_zerr
    err3 = [coeff["Berr"][2]]*n
    err3_s = pd.Series(data=err3)
    err_iz = df[key0["i_zerr"]]
    df[key1["i_zerr"]] = adderr_series(
      err2_s, err3_s, (coeff["B1"][1]-coeff["B1"][2])*err_gr, err_iz)

    return df
# Colors ======================================================================


# Reflectance =================================================================
def renormalize_ref(df, wnorm1, key_w, key_ref, key_referr):
    """Renormalize spectra at arbitrary wavelength (wnorm1).

    Parameters
    ----------
    df : pandas.DataFrame
        input data frame with reflectance
    wnorm1 : float
        wavelength at which the spectrum is normalized
    key_w : str
        keyword for wavelength
    key_ref : str
        keyword for reflectance
    key_referr : str
        keyword for reflectance error

    Return
    ------
    df : pandas.DataFrame
        output data frame with renormalized reflectance
    """
    
    # Check the wavelength is in the key_w
    w_list = df[key_w].values.tolist()
    if wnorm1 in w_list:
        print(f"Normalization wavelength {wnorm1} exists in the df.")
        pass
    else:
        print(f"  Normalization wavelength {wnorm1} does not exist in the df.")
        # Search the closest wavelength
        df["wdiff"] = abs(df[key_w] - wnorm1)
        idx_norm1 = df["wdiff"].idxmin()
        wnorm1 = df.at[idx_norm1, key_w] 
        print(f"  Normalize with the closest wavelength {wnorm1:.3f}")

    # Index of new normalized wavelength (wnorm1)
    idx_norm1 = df[df[key_w] == wnorm1].index.values[0]
    # Index of old norm where reflectance equals 1
    try:
        idx_norm0 = df[df[key_ref] == 1].index.values[0]
    except:
        df["refdiff"] = abs(df[key_ref] - 1)
        idx_norm0 = df["refdiff"].idxmin()

    print("")
    print(f"  Original norm. (w, idx) = ({df.at[idx_norm0, key_w]:.3f}, {idx_norm0})")
    print(f"  New      norm. (w, idx) = ({wnorm1:.3f}, {idx_norm1})")
    

    # Normalize reflectacne with reflectance at new wavelength
    df[key_referr] = df[key_referr]/df.at[idx_norm1, key_ref]
    df[key_ref] = df[key_ref]/df.at[idx_norm1, key_ref]
    return df


def calc_ref_from_color(df_in, sys, center="g"):
    """Calculate color reflectance from color.

    See DeMeo et al.2013 for detail description.  
    Typical solar color formats are below.

    1. Ivezic+2001 
        (u-g, g-r, r-i, i-z) = (1.32, 0.45, 0.10, 0.04)
    2. Holmberg+2006, DeMeo+2013
        (g-r, r-i, i-z)      = (0.45+-0.02, 0.12+-0.01, 0.04+-0.02)
        (g-g, r-g, i-g, z-g) = (0+-0, -0.45+-0.02, -0.55+-0.03, -0.61+-0.04)
    Both colors are in SDSS system.
    3. Willmer 2018 (We use this)
        SDSS (u, g, r, i, z) = (6.39, 5.11, 4.65, 4.53, 4.50)
        PS   (g, r, i, z, Y) = (5.03, 4.64, 4.52, 4.51, 4.50)

    Parameters
    ----------
    df_in : pandas.DataFrame
        dataframe contains magnitude/color
    sys : str
        magnitude system (SDSS or PS)
    center : str
        center band

    Return
    ------
    df : pandas.DataFrame
        dataframe contains color reflectance
        ref_g, referr_g, etc. are added.
    """
    df = df_in.copy()

    # Consistency =============================================================
    # 1. Choi+2023 refers to Willner (2018) and colors of the Sun in the 
    #    Pan-STARRS are (g−r) =0.61, (r−i) =0.35, and (i−z) =0.16.
    #    This is consistent with the SDSS & 'Vega' system in Willner 2018.
    #      g-r, r-i, i-z = 0.610, 0.350, 0.160
    #    (The selection of the system is wrong. Also, "Vega" system is wrong...) 

    # 2. Ye+2019 refers to Willner (2018) and colors of the Sun in the 
    #    Pan-STARRS are (g−r)=0.39.
    #    This is consistent with the Pan-STARRS & 'AB' system in Willner 2018.
    #      g-r = 0.39

    # 3. DeMeo & Carry (2013) refers to Holmberg et al. (2006) and colors of 
    #    the Sun in the SDSS are (r-g) = -0.45, (i-g) = -0.55, (z-g) = -0.61.
    #    This is consistent with the SDSS & "AB" system in Willmer 2018 
    #    (r-g, i-g, z-g = -0.460, -0.580, -0.610.)
    #    Bus it slightly differs.
    #
    # 4. Bolin et al. (2023) refers to Willner (2018) and colors of the Sun in the 
    #    SDSS are (g−r)=0.46 and (g-i)=0.58.
    #    This is consistent with the SDSS & "AB" system in Willmer 2018 
    # Consistency =============================================================

    # -> I confirmed that AB is correct, not Vega.
    
    # Note: g- and r-band normalization are available
    magsys = "AB"
    u_SDSS  = sun_mag_W18("u_SDSS", magsys)
    g_SDSS  = sun_mag_W18("g_SDSS", magsys)
    r_SDSS  = sun_mag_W18("r_SDSS", magsys)
    i_SDSS  = sun_mag_W18("i_SDSS", magsys)
    z_SDSS  = sun_mag_W18("z_SDSS", magsys)
    g_PS  = sun_mag_W18("g_PS", magsys)
    r_PS  = sun_mag_W18("r_PS", magsys)
    i_PS  = sun_mag_W18("i_PS", magsys)
    z_PS  = sun_mag_W18("z_PS", magsys)
    Y_PS  = sun_mag_W18("Y_PS", magsys)

    # Uncertainties are 0.02 to 0.04.
    magerr = 0.02
    uerr_SDSS, gerr_SDSS, rerr_SDSS, ierr_SDSS, zerr_SDSS = (
        magerr, magerr, magerr, magerr, magerr)
    gerr_PS, rerr_PS, ierr_PS, zerr_PS, Yerr_PS = (
        magerr, magerr, magerr, magerr, magerr)
    
    print("AB magnitude")
    print("Colors of the Sun in the Pan-STARRS system:")
    print(f"g-r, r-i, i-z = {g_PS - r_PS:.3f}, {r_PS - i_PS:.3f}, {i_PS - z_PS:.3f}")
    print("Colors of the Sun in the SDSS system:")
    print(f"r-g, i-g, z-g = {r_SDSS - g_SDSS:.3f}, {i_SDSS - g_SDSS:.3f}, {z_SDSS - g_SDSS:.3f}")

    # Solar color in PS
    if sys == "PS":
        if center == "g":
            g_g, g_gerr = 0, 0
            r_g, r_gerr = r_PS-g_PS, adderr(rerr_PS, gerr_PS)
            i_g, i_gerr = i_PS-g_PS, adderr(ierr_PS, gerr_PS)
            z_g, z_gerr = z_PS-g_PS, adderr(zerr_PS, gerr_PS)
            Y_g, Y_gerr = Y_PS-g_PS, adderr(Yerr_PS, gerr_PS)

            colors = ["g_g", "r_g", "i_g", "z_g", "Y_g"]
            c_sun = dict(
                g_g=0, 
                r_g=r_g, 
                i_g=i_g, 
                z_g=z_g, 
                Y_g=Y_g,
                g_gerr=0, 
                r_gerr=r_gerr, 
                i_gerr=i_gerr, 
                z_gerr=z_gerr, 
                Y_gerr=Y_gerr,
                )

        elif center == "r":
            g_r, g_rerr = g_PS-r_PS, adderr(gerr_PS, rerr_PS)
            r_r, r_rerr = 0, 0
            i_r, i_rerr = i_PS-r_PS, adderr(ierr_PS, rerr_PS)
            z_r, z_rerr = z_PS-r_PS, adderr(zerr_PS, rerr_PS)
            Y_r, Y_rerr = Y_PS-r_PS, adderr(Yerr_PS, rerr_PS)

            colors = ["g_r", "r_r", "i_r", "z_r", "Y_r"]
            c_sun = dict(
                g_r=g_r, 
                r_r=0, 
                i_r=i_r, 
                z_r=z_r, 
                Y_r=Y_r,
                g_rerr=g_rerr, 
                r_rerr=0, 
                i_rerr=i_rerr, 
                z_rerr=z_rerr, 
                Y_rerr=Y_rerr,
                )
        else:
            assert False, "Not implemented."
    else:
        assert False, "Not implemented."

    # Use common colors
    cols_df = df.columns.tolist()
    colors = list(set(cols_df) & set(colors))

    # Calculate reflectance
    for c in colors:
        x = -0.4*(df[c]-(c_sun[c]))
        print(f"  {c} (Sun, df) = ({c_sun[c]:.3f}, {float(df[c]):.3f})")
        xerr = 0.4*adderr(float(df[f"{c}err"]), c_sun[f"{c}err"])
        ref = 10**(x)
        # x can be negative (for blue object)
        referr = 10**x*np.log(10)*xerr
        # g, r, i or z
        band = c[0]
        df[f"ref_{band}"] = ref
        df[f"referr_{band}"] = referr
        print(f"  {band} (ref, referr) = ({float(ref):.3f}, {float(referr):.3f})")
        print("")

    return df


def renormalize_ref_Mahlke(
    df, wnorm1, key_w, key_ref, key_referru, key_referrl):
    """Renormalize spectra at arbitrary wavelength (wnorm1).

    Parameters
    ----------
    df : pandas.DataFrame
        input data frame with reflectance
    wnorm1 : float
        wavelength at which the spectrum is normalized
    key_w : str
        keyword for wavelength
    key_ref : str
        keyword for reflectance
    key_referru : str
        keyword for upper reflectance error 
    key_referrl : str
        keyword for lower reflectance error 

    Return
    ------
    df : pandas.DataFrame
        output data frame with renormalized reflectance
    """
    
    # Check the wavelength is in the key_w
    w_list = df[key_w].values.tolist()
    if wnorm1 in w_list:
        print(f"Normalization wavelength {wnorm1} exists in the df.")
        pass
    else:
        print(f"  Normalization wavelength {wnorm1} does not exist in the df.")
        # Search the closest wavelength
        df["wdiff"] = abs(df[key_w] - wnorm1)
        idx_norm1 = df["wdiff"].idxmin()
        wnorm1 = df.at[idx_norm1, key_w] 
        print(f"  Normalize with the closest wavelength {wnorm1:.3f}")

    # Index of new normalized wavelength (wnorm1)
    idx_norm1 = df[df[key_w] == wnorm1].index.values[0]
    # Index of old norm where reflectance equals 1
    try:
        idx_norm0 = df[df[key_ref] == 1].index.values[0]
    except:
        df["refdiff"] = abs(df[key_ref] - 1)
        idx_norm0 = df["refdiff"].idxmin()

    print(f"  Original norm. (w, idx) = ({df.at[idx_norm0, key_w]:.3f}, {idx_norm0})")
    print(f"  New      norm. (w, idx) = ({wnorm1:.3f}, {idx_norm1})")
    print("")
    

    # Normalize reflectacne with reflectance at new wavelength
    df[key_ref] = df[key_ref]/df.at[idx_norm1, key_ref]
    
    df[key_referru] = df[key_referru]/df.at[idx_norm1, key_ref]
    df[key_referrl] = df[key_referrl]/df.at[idx_norm1, key_ref]
    return df
# Reflectance =================================================================

def diff_nsigma(x1, err1, x2, err2):
    """Return difference in n-sigma
    """
    nsigma = np.abs(x1 - x2) / np.sqrt(err1**2 + err2**2)
    return nsigma
