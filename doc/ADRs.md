# Updating WCS in Headers 11/14/23

## Context

In updating WCS information there was a worry that information about datasets used in research would be lost.

## Decision

In order to ensure reproducibility and the ability to share consistent datasets "headerlets" were created. These hold the wcs inforamtion without any data.  

## Consequences

Headerlets can now be used to easily update and select WCS for datasets. Headerlets are files with the _hlet.fits suffix


# Updating APERTURE keyword through poller 11/14/23

## Context

Two types of poller files exist. Full poller files, that include all of the necessary information, and simple poller files that include only the filenames; additional information is taken from the header. We needed a way to update a header keyword ("APERTURE") and it was found that the easiest way to do this was to pass the keyword through the second column of the poller file.

## Decision

The second column of a poller file is now reserved (for WFPC2) for passing the aperture keyword through the code to update the header keyword. For other instruments, the code treats a two-column poller file as a simple (filename-only) poller file. 

## Consequences

Caution must be taken is using variations of the poller file while processing WFPC2 data.

# The Use of the TEAL Interface 11/14/23

## Context

The code can be run using the interactive GUI TEAL, however, it is hard to maintain. TEAL is also used in the code to load the json parameter files. 

## Decision

In order to make the code more easily maintainable, we will no longer support the use of the GUI for running tasks, however, TEAL will still be used in the background to load the parameter files as there is no current alternative for parsing that data. 

## Consequences

TEAL will need to be included in drizzlepac until a replacement for parsing the json files can be found. 