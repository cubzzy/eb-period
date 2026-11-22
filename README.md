**Code for period-finding for eclipsing binaries of specific spectral classes, along with random spectra obtained from SIMBAD**

We gathered the RA and Decs of eclipsing binaries using SIMBAD, obtained the light curves using ZTF's Multi-Object Search, then use this script to automatically compute the periods using Lomb-Scargle. 

**ipac_converter.py**

Takes the output from SIMBAD TAP query and converts it to an IPAC table uploadable to ZTF

**ztf_period.py**

Takes all the light curve data and automatically computes the periods, 2 times the period, and generates plots for Lomb-Scargle-1 and 2. Green highlighted plots are ones the script determined to be the true period (We found Lomb-Scargle-2 to be more accurate), but we still went through and selected ones that looked best to us, as it can be difficult to differentiate the primary and secondary periods with algorithms.

**make_histogram.py**

Makes a histogram of the periods based on spectral class and a histogram of the randomly selected binaries' periods. 

**Note** 

The plots show the pipeline's chosen period, but we hand-picked six periods (B Data6, 10, 15, 18, and 25, plus G Data29) which differed from the plots by 2 times the period because Lomb-Scargle had trouble differentiating the primary and secondary eclipses. 
Several of the scripts used spreadsheets that we manually built and are not included in this repository; this work was made for a class project and if you'd like to get more information or get the csv's we used reach out to supercalifornia1992@proton.me.