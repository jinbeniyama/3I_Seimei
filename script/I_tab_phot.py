#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Create a table about observing condition.
The result of query to JPL is saved just in case.
Output physical properties using JPL data as well.
"""
import os
from argparse import ArgumentParser as ap
import datetime
import pandas as pd
import numpy as np
from astropy.time import Time
from astroquery.jplhorizons import Horizons

from I_common import loc_Seimei, add_obsinfo, tab_phot


if __name__ == "__main__":
    parser = ap(description='Obtain 3I/ATLAS info. and create table.')
    parser.add_argument(
        "res", type=str, 
        help="Photometric results with")
    parser.add_argument(
        "--date", type=str, default=None,
        help="date to reuse the aspect data")
    parser.add_argument(
        "--outdir", type=str, default="tab",
        help="output directory")
    args = parser.parse_args()

    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    today = datetime.date.today()
    df_res = pd.read_csv(args.res, sep=" ")
    frameid_list = sorted(list(set(df_res.frameid)))
    print("Frame ids (original):")
    print(f"  {frameid_list}")
    assert len(frameid_list) == 12, "Check input file."
    # 12 sets
    # Exposure ID:
    #     gri 58521, 58528, 58534, 58541, 58542, 58543, 58544
    #     grz 58524, 58530, 58536, 58538, 58539
    # Actual orders:
    #  1. gri 
    #  2. grz star overlap
    #  3. gri
    #  4. grz
    #  5. gri
    #  6. grz 
    #  7. grz 
    #  8. grz
    #  9. gri 
    # 10. gri 
    # 11. gri 
    # 12. gri star overlap
    frameid_list = [
        "gri1", "grz1", "gri2", "grz2", "gri3", "grz3", 
        "grz4", "grz5", "gri4", "gri5", "gri6", "gri7"]
    print("Frame ids (sorted by hand):")
    print(f"  {frameid_list}")

    # Extract one arbitrary radius
    rad = np.min(df_res["radius"])
    df_res = df_res[df_res["radius"] == rad]

    if not args.date:
        # Seimei results ======================================================

        band_list = []
        jd_list = []
        utc_list = []
        for idx_obs, frameid in enumerate(frameid_list):
            # Use r-band
            df_temp = df_res[
                (df_res["frameid"] == frameid) & (df_res["band"] == "r")]
            df_temp = df_temp.reset_index(drop=True)
            assert len(df_temp) ==1, "Check the code." 

            # Starting time
            jd0 = np.min(df_temp["jd"])
            # Convert to utc
            utc0 = Time(
                str(jd0), format='jd', scale='utc').datetime.strftime('%Y-%m-%d %H:%M:%S')

            bands = frameid[0:3]
            band_list.append(bands)
            utc_list.append(utc0)
            jd_list.append(jd0)

        # r-band, measured by J.B. (Seimei)
        note_list = [
            "", "Overlap with a source", "", "", "",
            "", "", "", "", "",
            "", "Overlap with a source", 
            ]

        # Create dataframe
        df = pd.DataFrame(dict(
            obsdate0=utc_list,
            jd0=jd_list,
            band=band_list, 
            note=note_list,
            ))
        df["t_exp"] = 600
        # JPL style
        df["obj"] = "3I"
        
        # Add obs info.
        df = add_obsinfo(df, loc_Seimei)
        df = df.reset_index(drop=True)
        # Seimei results ======================================================

        out = f"I_phot_with_JPL_{today}.txt"
        out = os.path.join(args.outdir, out)
        df.to_csv(out, sep=" ", index=False)
        date = str(today)
    else:
        print(f"Reuse aspect data obtained on {args.date}.")
        f = f"I_phot_with_JPL_{args.date}.txt"
        f = os.path.join(args.outdir, f)
        df = pd.read_csv(f, sep=" ")
        date = args.date

    # create obs tex table
    out = f"I_tab_phot.tex"
    out = os.path.join(args.outdir, out)
    tab_phot(df, date, out)

    # Calculate distance in km
    obj = "3I"
    def angular_dist(distance_au, dist_arcsec):
        # 1 au in km
        AU_KM = 149597870.7  
        # 1 arcsec in radians
        arcsec_to_rad = np.deg2rad(1 / 3600) 
        dist_km = distance_au * AU_KM
        dist_rad = dist_arcsec * arcsec_to_rad
        size_km = dist_km * dist_rad
        return size_km

    ## Seimei
    t = "2025-07-15T14:00:00"
    jd = Time(str(t), format='isot', scale='utc').jd
    loc = loc_Seimei
    jpl = Horizons(id=obj, location=loc, epochs=jd)
    eph = jpl.ephemerides()
    delta = eph["delta"][0]
    r_h = eph["r"][0]
    ## In arcsec
    rad = 2.1
    dist_km = angular_dist(delta, rad)
    print(f"Heliocentric distance {r_h:.2f} au")
    print(f"Angular distance {dist_km:.1f} km")
    print()

    ## Bolin+2025
    t = "2025-07-06T15:00:00"
    jd = Time(str(t), format='isot', scale='utc').jd
    loc = 500
    jpl = Horizons(id=obj, location=loc, epochs=jd)
    eph = jpl.ephemerides()
    delta = eph["delta"][0]
    r_h = eph["r"][0]
    ## In arcsec
    rad = 4
    dist_km = angular_dist(delta, rad)
    print(f"Heliocentric distance {r_h:.2f} au")
    print(f"Angular distance {dist_km:.1f} km")
