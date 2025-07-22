# Spectrophotometry of 3I/ATLAS
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[developer mail](mailto:jbeniyama@oca.eu)

## Overview
This is a repository for 3I/ATLAS paper in prep.
Figures are made in `./plot`.

## Analysis (in /analysis)
Do photometry and derive magnitudes and colors of 3I. The results (`mag_3I_20250715.txt` and `col_3I_20250715.txt`) are saved in `/data`
* `photometry_3I.ipynb`: Do photometry of 3I (and Popescu for validation)
* `derive_mag_col_3I.ipynb`: Derive magnitudes and colors of 3I (and Popescu for validation)

## Data (in /data)
* `mag_3I_20250715.txt`: Calibrated magnitude of 3I/ATLAS on 2025-07-15
* `col_3I_20250715.txt`: Calibrated colors of 3I/ATLAS on 2025-07-15

## Plotting figures in the paper (hit commands in /plot)
```
# Obtain ephem
python3 ../script/I_fig_ephem.py --obtain
# Just plot
python3 ../script/I_fig_ephem.py
```

```
# Table 
python3 ../script/I_tab_phot.py ../../data/photres_3I_mean_600s_20250715.txt
```

```
# Cutout
python3 ../script/I_fig_cutout.py
```

```
# Lightcurvea
python3 ../script/I_fig_lc.py ../data/3I_mag_20250715.txt
```

```
# Color-color diagram
python3 ../script/I_fig_cc.py ../data/col_3I_20250715.txt ../data/col_Popescu_20250715.txt 6 --outtype pdf --magtype PS
```

```
# Reflectance
python3 ../script/I_fig_ref.py ../data/col_3I_20250715.txt 6 --outtype pdf
```

## Dependencies
This repository is depending on `Python`, `NumPy`, `SciPy`, `SEP`, `Astropy`, `Astroquery`.
Scripts are developed on `Python 3.9.15`, `NumPy 1.25.0`, `SciPy 1.11.1`, `SEP 1.2.0`, `Astropy 5.3`, `Astroquery 0.4.6`.
