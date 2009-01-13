#!/usr/bin/env python
"""

Function for computing and subtracting the backgroud of an image.  The
algorithm employed here uses a sigma clipped median of  each *sci* image in a data file.  
Then the sky value for each detector
is compared and the lowest value is subtracted from all chips in the
detector.  Finally, the MDRIZSKY keyword is updated in the header
of the input files.

:author: Christopher Hanley
:author: Megan Sosey

"""

def subtract(inputlist, 
             skywidth,
             skystat,
             skylower,
             skyupper,
             skyclip,
             skylsignma,
             skyusignma,
             skyuser
             ):

    # Print out the parameters provided by the interface
    print "USER PARAMETERS:"
    print "skywidth  = ",skywidth
    print "skystat   = ",skystat
    print "skylower  = ",skylower
    print "skyupper  = ",skyupper
    print "skyclip   = ",skyclip
    print "skylsigma = ",skylsigma
    print "skyusigma = ",skyusigma
    print "skyuser   = ",skyuser
    print "\n"

    """ Processes sky in input images."""
    if skyuser != '':
        # User Subtraction Case, User has done own sky subtraction, we use 
        # image header value for subtractedsky value
        print "User has done own sky subtraction, updating MDRIZSKY with supplied value..."
        for p in self.assoc.parlist:
            if int(p['group']) == 1:
                try:
                    handle = fileutil.openImage(p['image'].datafile,mode='update',memmap=0)
                    userSkyValue = handle[0].header[skypars['skyuser']]
                    handle.close()
                except:
                    print "**************************************************************"
                    print "*"
                    print "*  Cannot find keyword ",skypars['skyuser']," in ",p['image'].datafile," to update"
                    print "*"
                    print "**************************************************************\n\n\n"
                    raise KeyError
            p['image'].setSubtractedSky(userSkyValue)
    elif (skysub):
        # Compute our own sky values and subtract them from the image copies.
        print "Subtracting sky..."
        imageMinDict = {}
        currentImageName = "no match"
        for p in self.assoc.parlist:
            p['image'].computeSky(skypars)
            computedImageSky = p['image'].getComputedSky()
            if (p['rootname'] != currentImageName):
                currentMinSky = computedImageSky
                currentImageName = p['rootname']
                imageMinDict[currentImageName]=currentMinSky
            else:
                if (computedImageSky < imageMinDict[p['rootname']]):
                    imageMinDict[p['rootname']]=computedImageSky

        for p in self.assoc.parlist:
            # We need to account for the fact that STIS associations contain
            # separate exposures for the same chip within the same file.
            # In those cases, we do not want to use the minimum sky value
            # for the set of chips, but rather use each chip's own value.
            # NOTE: This can be generalized later with changes to PyDrizzle
            #       to provide an attribute that specifies whether each member
            #       associated with file is a separate exposure or not.
            #   WJH/CJH 
            if (p['exposure'].header['INSTRUME'] != 'STIS'):
                p['image'].setSubtractedSky(imageMinDict[p['rootname']])
            else:
                p['image'].setSubtractedSky(p['image'].getComputedSky())
            p['image'].subtractSky()

    else:
        # Default Case, sky subtraction is turned off.  No sky subtraction done to image.
        print "No sky subtraction requested, MDRIZSKY set to a value of 0..."
        for p in self.assoc.parlist:
            p['image'].setSubtractedSky(0)

    #Update the MDRIZSKY KEYWORD
    for p in self.assoc.parlist:
        if int(p['group']) == 1:
            p['image'].updateMDRIZSKY(p['orig_filename'])
            #p['image'].updateMDRIZSKY(p['image'].datafile)
