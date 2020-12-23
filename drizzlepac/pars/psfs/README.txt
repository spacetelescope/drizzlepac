The PSFs in this directory all have filenames based on:

    <instrument>_<detector>_<filter>_psf.fits

where all terms are lower-case.  In addition, each PSF will represent a single
filter, not a combination of filters.  

Each PSF was created using these parameters:
  * chip 2 (for WFC3/UVIS, ACS/WFC)
  * position: center of full chip ( [512, 512] or [2048, 1024] )
  * Spectra:  Blackbody (for wide- and medium-band filters)
    - T=10000K for filters < 4000A 
    - T=5000K for filters < 4000A
  * secondary mirror positions: +/- 1.0 microns
