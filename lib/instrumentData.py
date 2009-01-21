#!/usr/bin/env python
"""
functions associated with HST instruments


"""

import sys
import pytools
import util
import acsdata


def getBasicInstrData(filename):
    """ return a dictionary with basic instrument data taken from the header
   		of the input filename
    """
        
    imageHandle=fileutil.openImage(filename,mode='update',memmap=0)
    priHeader=imageHandle[0].header
	
    instrument=priHeader["INSTRUME"]
    
    #these fill in specific instrument information
    #instrData is returned as a dictionary of key vals
    if ("ACS" in instrument):
        instrData=getACSInfo(priHeader)
    if ("NICMOS" in instrument):
        instrData=getNICMOSInfo(priHeader)
    if ("WFC3" in instrument):
        instrData=getWFC3Info(priHeader)
    if("WFPC2" in instrument):
        instrData=getWFPC2Info(priHeader)

        
    #keywords which are common to all instruments
    genericKW=["INSTRUME","NAXIS1","NAXIS2","LTV1","LTV2","EXPTIME","NEXTEND"]
    
    for key in genericKW:
    	instrData[key]=priHeader[key]

    
    imageHandle.close()
        
    return instrData    
    
    
class InputImage:
    '''The InputImage class is the base class for all of the various
       types of images
    '''

    def __init__(self, filename):
        self.filename = filename
        self.rootname = util.findrootname(filename)
        self.subtractedSky=0.0 #sky subtracted from all the chips for the instrument
        
        setInstrumentParameters(self)
        
        
        
    def setInstrumentParameters(self, instrpars, pri_header):
        """ 
        Sets the instrument parameters.
        """
        self.refplatescale=0.0
        self.instrumentName=None
        self.numberOfChips=1 #these are directly related to the number of science images in the file
        self.dataUnits="electrons" #set to the units the science data is in, as read from header
        
        pass
    
    def doUnitConversions(self):
        """
        Convert the sci extension pixels to electrons
        """
        pass
        
    def getInstrParameter(self, value, header, keyword):
        """ This method gets a instrument parameter from a
            pair of task parameters: a value, and a header keyword.

            The default behavior is:
              - if the value and header keyword are given, raise an exception.
              - if the value is given, use it.
              - if the value is blank and the header keyword is given, use
                the header keyword.
              - if both are blank, or if the header keyword is not
                found, return None.
        """
        if (value != None and value != '')  and (keyword != None and keyword.strip() != ''):
            exceptionMessage = "ERROR: Your input is ambiguous!  Please specify either a value or a keyword.\n  You specifed both " + str(value) + " and " + str(keyword) 
            raise ValueError, exceptionMessage
        elif value != None and value != '':
            return self._averageFromList(value)
        elif keyword != None and keyword.strip() != '':
            return self._averageFromHeader(header, keyword)
        else:
            return None

    def _averageFromHeader(self, header, keyword):
        """ Averages out values taken from header. The keywords where
            to read values from are passed as a comma-separated list.
        """
        _list = ''
        for _kw in keyword.split(','):
            if header.has_key(_kw):
                _list = _list + ',' + str(header[_kw])
            else:
                return None
        return self._averageFromList(_list)

    def _averageFromList(self, param):
        """ Averages out values passed as a comma-separated
            list, disregarding the zero-valued entries.
        """
        _result = 0.0
        _count = 0

        for _param in param.split(','):
            if _param != '' and float(_param) != 0.0:
                _result = _result + float(_param)
                _count  += 1

        if _count >= 1:
            _result = _result / _count
        return _result

    def getreferencesky(self):
        return (self._subtractedsky * (self.refplatescale / self.platescale)**2 )                
                
    def updateMDRIZSKY(self,filename=None):
    
        if (filename == None):
            filename = self.name
            
        try:
            _handle = fileutil.openImage(filename,mode='update',memmap=self.memmap)
        except:
            raise IOError, "Unable to open %s for sky level computation"%filename
        try:
            try:
                # Assume MDRIZSKY lives in primary header
                print "Updating MDRIZSKY in %s with %f"%(filename,self.getSubtractedSky())
                _handle[0].header['MDRIZSKY'] = self.getSubtractedSky()
            except:
                print "Cannot find keyword MDRIZSKY in %s to update"%filename
                print "Adding MDRIZSKY keyword to primary header with value %f"%self.getSubtractedSky()
                _handle[0].header.update('MDRIZSKY',self.getSubtractedSky(), 
                    comment="Sky value subtracted by Multidrizzle")
        finally:
            _handle.close()
