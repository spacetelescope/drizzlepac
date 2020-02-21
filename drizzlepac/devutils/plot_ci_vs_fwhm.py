"""
1: read in ci_fwhm file(s)
2: smash all input datasets into a single dataset
3: determine min, max CI values
4: set up CI value bins based on dataset min/max values and specified bin size
5: for each bin, compute resistant mean FWHM value, sigma value
6: plot bin CI value vs. mean FWHM values
"""