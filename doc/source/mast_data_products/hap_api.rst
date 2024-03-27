Hubble Advanced Products API
============================
.. _hap_api:

Drizzlepac has evolved over time from not only being responsible for creating the calibration pipeline
drizzled products (i.e., ipppssoot_dr[c|z].fits), but also to being responsible for the generation of
the two types of Hubble Advanced Products (HAP). The HAP products are Single Visit Mosaics (SVMs) and
Multiple Visit Mosaics (MVMs) and the creation of these products has been incorporated into the
calibration pipeline processing. SVM and MVM drizzled products are distingushed from pipeline
drizzled counterparts and each other through their filenames, as well as internal contents.

SVM: hst_<propid>_<obsetid>_<instr>_<detector>_<filter>_<ipppssoo>_dr[c|z].fits)
MVM: hst_skycell_p<PPPP>x<XX>y<YY>_<instr>_<detector>_<filter>_<label>_dr[c|z].fits)

where PPPP = 4-digit projection cell number
XX, YY = 2-digit sky cell numbers

runastrodriz.py is the module to control operation of AstroDrizzle to remove distortion and combine HST images.

runsinglehap.py is the module which controls the SVM processing.

runmultihap.py is the module which controls the MVM processing.

Sucessful code will result in an exit code of 0.  If there is an error, the exit code will be non-zero. 
An exit code of 55 is for data containing only SBC or HRC data, and an exit code of 65 is for "no viable data".

.. toctree:: 
    :maxdepth: 1

    astrometry_api
    hap-parameters
    run_singlehap
    mvm_api
    run_multihap
    mvm_utilities_api
    make_custom_mosaic
    object_classes
    