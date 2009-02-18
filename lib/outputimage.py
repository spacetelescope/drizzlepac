import types
import pyfits
from pytools import fileutil, readgeis

yes = True
no = False

RESERVED_KEYS = ['NAXIS','BITPIX','DATE','IRAF-TLM','XTENSION','EXTNAME','EXTVER']

EXTLIST = ('SCI', 'WHT', 'CTX')

DTH_KEYWORDS=['CD1_1','CD1_2', 'CD2_1', 'CD2_2', 'CRPIX1',
'CRPIX2','CRVAL1', 'CRVAL2', 'CTYPE1', 'CTYPE2']

class OutputImage:
    """
        This class manages the creation of the array objects
        which will be used by Drizzle. The three arrays, SCI/WHT/CTX,
        will be setup either as extensions in a
        single multi-extension FITS file, or as separate FITS
        files.
    """
    def __init__(self, plist, input_pars, build=yes, wcs=None, single=no, blot=no):
        """
        The object 'plist' must contain at least the following members:
            plist['output']  - name of output FITS image (for SCI)
            plist['outnx']    - size of X axis for output array
            plist['outny']    - size of Y axis for output array
        If 'single=yes', then 'plist' also needs to contain:
            plist['outsingle']
            plist['outsweight']
            plist['outscontext']
        If 'blot=yes', then 'plist' also needs:
            plist['data']
            plist['outblot']
            plist['blotnx'],plist['blotny']

        If 'build' is set to 'no', then each extension/array must be
        in separate FITS objects.  This would also require:
          plist['outdata']    - name of output SCI FITS image
          plist['outweight']  - name of output WHT FITS image
          plist['outcontext'] - name of output CTX FITS image

        Optionally, the overall exposure time information can be passed as:
            plist['texptime'] - total exptime for output
            plist['expstart'] - start time of combined exposure
            plist['expend']   - end time of combined exposure


        """
        self.build = build
        self.single = single
        self.parlist = plist
        self.input_pars = input_pars
        _nimgs = len(self.parlist)
        self.bunit = None
        self.units = 'cps'
        
        # Merge input_pars with each chip's outputNames object
        for p in self.parlist:
            p.update(input_pars)
            
        if not blot:
            self.output = plist[0]['output']
            self.shape = (plist[0]['outny'],plist[0]['outnx'])
        else:
            self.output = plist[0]['blotImage']
            self.shape = (plist[0]['blotny'],plist[0]['blotnx'])

        # Keep track of desired output WCS computed by PyDrizzle
        self.wcs = wcs

        #
        # Apply special operating switches:
        #   single - separate output for each input
        #
        if single:
            _outdata = plist[0]['outSingle']
            _outweight = plist[0]['outSWeight']
            _outcontext = plist[0]['outSContext']
            # Only report values appropriate for single exposure
            self.texptime = plist[0]['texptime']
            self.expstart = plist[0]['texpstart']
            self.expend = plist[0]['texpend']
        else:
            if self.build:
                _outdata = plist[0]['outFinal']
            else:
                _outdata = plist[0]['outSci']
            _outweight = plist[0]['outWeight']
            _outcontext = plist[0]['outContext']
            # Report values appropriate for entire combined product
            self.texptime = plist[0]['texptime']
            self.expstart = plist[0]['texpstart']
            self.expend = plist[_nimgs-1]['texpend']


        if blot:
            _outdata = plist[0]['blotImage']

        if not self.build or single:
            self.output = _outdata

        self.outdata = _outdata
        self.outweight = _outweight
        self.outcontext = _outcontext

    def set_bunit(self,bunit):
        """ Method used to update the value of the bunit attribute."""
        self.bunit = bunit
        
    def set_units(self,units):
        """ Method used to record what units were specified by the user
        for the output product."""
        self.units = units
        

    def writeFITS(self, template, sciarr, whtarr, ctxarr=None, versions=None, extlist=EXTLIST, overwrite=yes):
        """ Generate PyFITS objects for each output extension
            using the file given by 'template' for populating
            headers.

            The arrays will have the size specified by 'shape'.
        """

        if fileutil.findFile(self.output):
            if overwrite:
                print 'Deleting previous output product: ',self.output
                fileutil.removeFile(self.output)

            else:
                print 'WARNING:  Output file ',self.output,' already exists and overwrite not specified!'
                print 'Quitting... Please remove before resuming operations.'
                raise IOError

        # Default value for NEXTEND when 'build'== True
        nextend = 3
        if not self.build:
            nextend = 0
            if self.outweight:
                if overwrite:
                    if fileutil.findFile(self.outweight):
                        print 'Deleting previous output WHT product: ',self.outweight
                    fileutil.removeFile(self.outweight)
                else:
                    print 'WARNING:  Output file ',self.outweight,' already exists and overwrite not specified!'
                    print 'Quitting... Please remove before resuming operations.'
                    raise IOError


            if self.outcontext:
                if overwrite:
                    if fileutil.findFile(self.outcontext):
                        print 'Deleting previous output CTX product: ',self.outcontext
                    fileutil.removeFile(self.outcontext)
                else:
                    print 'WARNING:  Output file ',self.outcontext,' already exists and overwrite not specified!'
                    print 'Quitting... Please remove before resuming operations.'
                    raise IOError

        # Get default headers from multi-extension FITS file
        # If input data is not in MEF FITS format, it will return 'None'
        # and those headers will have to be generated from drizzle output
        # file FITS headers.
        # NOTE: These are HEADER objects, not HDUs
        prihdr,scihdr,errhdr,dqhdr = getTemplates(template,extlist)

        if prihdr == None:
            # Use readgeis to get header for use as Primary header.
            _indx = template.find('[')
            if _indx < 0:
                _data = template
            else:
                _data = template[:_indx]

            fpri = readgeis.readgeis(_data)
            prihdr = fpri[0].header.copy()
            fpri.close()
            del fpri


        # Setup primary header as an HDU ready for appending to output FITS file
        prihdu = pyfits.PrimaryHDU(header=prihdr,data=None)

        # Start by updating PRIMARY header keywords...
        prihdu.header.update('EXTEND',pyfits.TRUE,after='NAXIS')
        prihdu.header.update('NEXTEND',nextend)
        prihdu.header.update('FILENAME', self.output)

        # Update the ROOTNAME with the new value as well
        _indx = self.output.find('_drz')
        if _indx < 0:
            prihdu.header.update('ROOTNAME', self.output)
        else:
            prihdu.header.update('ROOTNAME', self.output[:_indx])


        # Get the total exposure time for the image
        # If not calculated by PyDrizzle and passed through
        # the pardict, then leave value from the template image.
        if self.texptime:
            prihdu.header.update('EXPTIME', self.texptime)
            prihdu.header.update('EXPSTART', self.expstart)
            prihdu.header.update('EXPEND', self.expend)

        #Update ASN_MTYPE to reflect the fact that this is a product
        # Currently hard-wired to always output 'PROD-DTH' as MTYPE
        prihdu.header.update('ASN_MTYP', 'PROD-DTH')

        # Update DITHCORR calibration keyword if present
        # Remove when we can modify FITS headers in place...
        if prihdu.header.has_key('DRIZCORR') > 0:
            prihdu.header['DRIZCORR'] = 'COMPLETE'
        if prihdu.header.has_key('DITHCORR') > 0:
            prihdu.header['DITHCORR'] = 'COMPLETE'


        prihdu.header.update('NDRIZIM',len(self.parlist),
            comment='Drizzle, No. images drizzled onto output')

        # Only a subset of these keywords makes sense for the new WCS based
        # transformations. They need to be reviewed to decide what to keep
        # and what to leave out.
        #self.addDrizKeywords(prihdu.header,versions)

        if scihdr:
            del scihdr['OBJECT']
            if scihdr.has_key('CCDCHIP'): scihdr.update('CCDCHIP','-999')
            if scihdr.has_key('NCOMBINE') > 0:
                scihdr.update('NCOMBINE', self.parlist[0]['nimages'])

            # If BUNIT keyword was found and reset, then 
            if self.bunit is not None:
                scihdr.update('BUNIT',self.bunit,comment="Units of science product")

            if self.wcs:
                # Update WCS Keywords based on PyDrizzle product's value
                # since 'drizzle' itself doesn't update that keyword.
                scihdr.update('ORIENTAT',self.wcs.orientat)
                scihdr.update('CD1_1',self.wcs.wcs.cd[0][0])
                scihdr.update('CD1_2',self.wcs.wcs.cd[0][1])
                scihdr.update('CD2_1',self.wcs.wcs.cd[1][0])
                scihdr.update('CD2_2',self.wcs.wcs.cd[1][1])
                scihdr.update('CRVAL1',self.wcs.wcs.crval[0])
                scihdr.update('CRVAL2',self.wcs.wcs.crval[1])
                scihdr.update('CRPIX1',self.wcs.wcs.crpix[0])
                scihdr.update('CRPIX2',self.wcs.wcs.crpix[1])
                scihdr.update('VAFACTOR',1.0)
                # Remove any reference to TDD correction
                if scihdr.has_key('TDDALPHA'):
                    del scihdr['TDDALPHA']
                    del scihdr['TDDBETA']
                # Remove '-SIP' from CTYPE for output product
                if scihdr['ctype1'].find('SIP') > -1:
                    scihdr.update('ctype1', scihdr['ctype1'][:-4])
                    scihdr.update('ctype2',scihdr['ctype2'][:-4])
                    # Remove SIP coefficients from DRZ product
                    for k in scihdr.items():
                        if (k[0][:2] in ['A_','B_']) or (k[0][:3] in ['IDC','SCD'] and k[0] != 'IDCTAB') or \
                        (k[0][:6] in ['SCTYPE','SCRVAL','SNAXIS','SCRPIX']): 
                            del scihdr[k[0]]

        ##########
        # Now, build the output file
        ##########
        if self.build:
            print '-Generating multi-extension output file: ',self.output
            fo = pyfits.HDUList()

            # Add primary header to output file...
            fo.append(prihdu)

            hdu = pyfits.ImageHDU(data=sciarr,header=scihdr,name=extlist[0])
            fo.append(hdu)

            # Build WHT extension here, if requested...
            if errhdr:
                errhdr.update('CCDCHIP','-999')

            hdu = pyfits.ImageHDU(data=whtarr,header=errhdr,name=extlist[1])
            hdu.header.update('EXTVER',1)
            if self.wcs:
                # Update WCS Keywords based on PyDrizzle product's value
                # since 'drizzle' itself doesn't update that keyword.
                hdu.header.update('ORIENTAT',self.wcs.orientat)
                hdu.header.update('CD1_1',self.wcs.wcs.cd[0][0])
                hdu.header.update('CD1_2',self.wcs.wcs.cd[0][1])
                hdu.header.update('CD2_1',self.wcs.wcs.cd[1][0])
                hdu.header.update('CD2_2',self.wcs.wcs.cd[1][1])
                hdu.header.update('CRVAL1',self.wcs.wcs.crval[0])
                hdu.header.update('CRVAL2',self.wcs.wcs.crval[1])
                hdu.header.update('CRPIX1',self.wcs.wcs.crpix[0])
                hdu.header.update('CRPIX2',self.wcs.wcs.crpix[1])
                hdu.header.update('VAFACTOR',1.0)

            fo.append(hdu)

            # Build CTX extension here
            # If there is only 1 plane, write it out as a 2-D extension
            if self.outcontext:
                if ctxarr.shape[0] == 1:
                    _ctxarr = ctxarr[0]
                else:
                    _ctxarr = ctxarr
            else:
                _ctxarr = None

            hdu = pyfits.ImageHDU(data=_ctxarr,header=dqhdr,name=extlist[2])
            hdu.header.update('EXTVER',1)
            if self.wcs:
                # Update WCS Keywords based on PyDrizzle product's value
                # since 'drizzle' itself doesn't update that keyword.
                hdu.header.update('ORIENTAT',self.wcs.orientat)
                hdu.header.update('CD1_1',self.wcs.wcs.cd[0][0])
                hdu.header.update('CD1_2',self.wcs.wcs.cd[0][1])
                hdu.header.update('CD2_1',self.wcs.wcs.cd[1][0])
                hdu.header.update('CD2_2',self.wcs.wcs.cd[1][1])
                hdu.header.update('CRVAL1',self.wcs.wcs.crval[0])
                hdu.header.update('CRVAL2',self.wcs.wcs.crval[1])
                hdu.header.update('CRPIX1',self.wcs.wcs.crpix[0])
                hdu.header.update('CRPIX2',self.wcs.wcs.crpix[1])
                hdu.header.update('VAFACTOR',1.0)
                

            fo.append(hdu)

            fo.writeto(self.output)
            fo.close()
            del fo, hdu

        else:
            print '-Generating simple FITS output: ',self.outdata
            fo = pyfits.HDUList()

            hdu = pyfits.PrimaryHDU(data=sciarr, header=prihdu.header)

            # Append remaining unique header keywords from template DQ
            # header to Primary header...
            if scihdr:
                for _card in scihdr.ascard:
                    if _card.key not in RESERVED_KEYS and hdu.header.has_key(_card.key) == 0:
                        hdu.header.ascard.append(_card)
            del hdu.header['PCOUNT']
            del hdu.header['GCOUNT']
            hdu.header.update('filename',self.outdata)

            # Add primary header to output file...
            fo.append(hdu)
            fo.writeto(self.outdata)
            del fo,hdu

            if self.outweight and whtarr != None:
                # We need to build new PyFITS objects for each WHT array
                fwht = pyfits.HDUList()

                if errhdr:
                    errhdr.update('CCDCHIP','-999')

                hdu = pyfits.PrimaryHDU(data=whtarr, header=prihdu.header)

                # Append remaining unique header keywords from template DQ
                # header to Primary header...
                if errhdr:
                    for _card in errhdr.ascard:
                        if _card.key not in RESERVED_KEYS and hdu.header.has_key(_card.key) == 0:
                            hdu.header.ascard.append(_card)
                hdu.header.update('filename',self.outweight)
                hdu.header.update('CCDCHIP','-999')
                if self.wcs:
                    # Update WCS Keywords based on PyDrizzle product's value
                    # since 'drizzle' itself doesn't update that keyword.
                    hdu.header.update('ORIENTAT',self.wcs.orientat)
                    hdu.header.update('CD1_1',self.wcs.wcs.cd[0][0])
                    hdu.header.update('CD1_2',self.wcs.wcs.cd[0][1])
                    hdu.header.update('CD2_1',self.wcs.wcs.cd[1][0])
                    hdu.header.update('CD2_2',self.wcs.wcs.cd[1][1])
                    hdu.header.update('CRVAL1',self.wcs.wcs.crval[0])
                    hdu.header.update('CRVAL2',self.wcs.wcs.crval[1])
                    hdu.header.update('CRPIX1',self.wcs.wcs.crpix[0])
                    hdu.header.update('CRPIX2',self.wcs.wcs.crpix[1])
                    hdu.header.update('VAFACTOR',1.0)

                # Add primary header to output file...
                fwht.append(hdu)
                fwht.writeto(self.outweight)
                del fwht,hdu

            # If a context image was specified, build a PyFITS object
            # for it as well...
            if self.outcontext and ctxarr != None:
                fctx = pyfits.HDUList()

                # If there is only 1 plane, write it out as a 2-D extension
                if ctxarr.shape[0] == 1:
                    _ctxarr = ctxarr[0]
                else:
                    _ctxarr = ctxarr

                hdu = pyfits.PrimaryHDU(data=_ctxarr, header=prihdu.header)

                # Append remaining unique header keywords from template DQ
                # header to Primary header...
                if dqhdr:
                    for _card in dqhdr.ascard:
                        if _card.key not in RESERVED_KEYS and hdu.header.has_key(_card.key) == 0:
                            hdu.header.ascard.append(_card)
                hdu.header.update('filename', self.outcontext)
                if self.wcs:
                    # Update WCS Keywords based on PyDrizzle product's value
                    # since 'drizzle' itself doesn't update that keyword.
                    hdu.header.update('ORIENTAT',self.wcs.orientat)
                    hdu.header.update('CD1_1',self.wcs.wcs.cd[0][0])
                    hdu.header.update('CD1_2',self.wcs.wcs.cd[0][1])
                    hdu.header.update('CD2_1',self.wcs.wcs.cd[1][0])
                    hdu.header.update('CD2_2',self.wcs.wcs.cd[1][1])
                    hdu.header.update('CRVAL1',self.wcs.wcs.crval[0])
                    hdu.header.update('CRVAL2',self.wcs.wcs.crval[1])
                    hdu.header.update('CRPIX1',self.wcs.wcs.crpix[0])
                    hdu.header.update('CRPIX2',self.wcs.wcs.crpix[1])
                    hdu.header.update('VAFACTOR',1.0)

                fctx.append(hdu)
                fctx.writeto(self.outcontext)
                del fctx,hdu


    def addDrizKeywords(self,hdr,versions):
        """ Add drizzle parameter keywords to header. """

        # Extract some global information for the keywords
        _geom = 'User parameters'

        _imgnum = 0
        for pl in self.parlist:

            # Start by building up the keyword prefix based
            # on the image number for the chip
            _imgnum += 1
            _keyprefix = 'D%03d'%_imgnum
            if not isinstance(pl['driz_mask'],types.StringType):
                _driz_mask_name = 'static mask'
            else:
                _driz_mask_name = pl['driz_mask']

            hdr.update(_keyprefix+'VER',pl['driz_version'][:44],
                comment='Drizzle, task version')

    #       Then the source of the geometric information
            hdr.update(_keyprefix+'GEOM','User parameters',
                comment= 'Drizzle, source of geometric information')

    #       Now we continue to add the other items using the same
    #       "stem"
            hdr.update(_keyprefix+'DATA',pl['data'][:64],
                comment= 'Drizzle, input data image')

            hdr.update(_keyprefix+'DEXP',pl['exptime'],
                comment= 'Drizzle, input image exposure time (s)')

            hdr.update(_keyprefix+'OUDA',pl['outFinal'][:64],
                comment= 'Drizzle, output data image')

            hdr.update(_keyprefix+'OUWE',pl['outWeight'][:64],
                comment= 'Drizzle, output weighting image')

            hdr.update(_keyprefix+'OUCO',pl['outContext'][:64],
                comment= 'Drizzle, output context image')

            hdr.update(_keyprefix+'MASK',_driz_mask_name[:64],
                comment= 'Drizzle, input weighting image')

            # Process the values of WT_SCL to be consistent with
            # what IRAF Drizzle would output
            if pl['wt_scl'] == 'exptime': _wtscl = pl['exptime']
            elif pl['wt_scl'] == 'expsq': _wtscl = pl['exptime']*pl['exptime']
            else: _wtscl = pl['wt_scl']

            hdr.update(_keyprefix+'WTSC',_wtscl,
                comment= 'Drizzle, weighting factor for input image')

            hdr.update(_keyprefix+'KERN',pl['kernel'],
                comment= 'Drizzle, form of weight distribution kernel')

            hdr.update(_keyprefix+'PIXF',pl['pixfrac'],
                comment= 'Drizzle, linear size of drop')

            hdr.update(_keyprefix+'COEF',"SIP",
                comment= 'Drizzle, coefficients file name ')

            hdr.update(_keyprefix+'XGIM',"SIP",
                comment= 'Drizzle, X distortion image name ')

            hdr.update(_keyprefix+'YGIM',"SIP",
                comment= 'Drizzle, Y distortion image name ')

            #hdr.update(_keyprefix+'LAM',pl['plam'],
            #    comment='Drizzle, wavelength applied for transformation (nm)')
            """
    #       Only put the next entries is we are NOT using WCS
            hdr.update(_keyprefix+'SCAL',pl['scale'],
             comment=   'Drizzle, scale (pixel size) of output image')

    #       Convert the rotation angle back to degrees
            hdr.update(_keyprefix+'ROT',float("%0.8f"%pl['rot']),
             comment= 'Drizzle, rotation angle, degrees anticlockwise')

    #       Check the SCALE and units
    #       The units are INPUT pixels on entry to this routine
            hdr.update(_keyprefix+'XSH',pl['xsh'],
             comment= 'Drizzle, X shift applied')

            hdr.update(_keyprefix+'YSH',pl['ysh'],
             comment= 'Drizzle, Y shift applied')

            hdr.update(_keyprefix+'SFTU','output',
             comment='Drizzle, units used for shifts')

            hdr.update(_keyprefix+'SFTF','output',
             comment=   'Drizzle, frame in which shifts were applied')

            hdr.update(_keyprefix+'EXKY','EXPTIME',
                comment= 'Drizzle, exposure keyword name in input image')

            hdr.update(_keyprefix+'INUN','counts',
                comment= 'Drizzle, units of input image - counts or cps')
            """
            hdr.update(_keyprefix+'OUUN',self.units,
                comment= 'Drizzle, units of output image - counts or cps')

            hdr.update(_keyprefix+'FVAL',pl['fillval'],
                comment= 'Drizzle, fill value for zero weight output pix')

            OFF=0.5

            hdr.update(_keyprefix+'INXC',float(pl['blotnx']/2)+OFF,
                comment= 'Drizzle, reference center of input image (X)')

            hdr.update(_keyprefix+'INYC',float(pl['blotny']/2)+OFF,
                comment= 'Drizzle, reference center of input image (Y)')

            hdr.update(_keyprefix+'OUXC',float(pl['outnx']/2)+OFF,
                comment= 'Drizzle, reference center of output image (X)')

            hdr.update(_keyprefix+'OUYC',float(pl['outny']/2)+OFF,
                comment= 'Drizzle, reference center of output image (Y)')

        # Add version information as HISTORY cards to the header
        if versions != None:
            ver_str = "PyDrizzle processing performed using: "
            hdr.add_history(ver_str)
            for k in versions.keys():
                ver_str = '    '+str(k)+' Version '+str(versions[k])
                hdr.add_history(ver_str)



def getTemplates(fname,extlist):
    # Obtain default headers for output file
    # If the output file already exists, use it
    # If not, use an input file for this information.
    #
    # NOTE: Returns 'pyfits.Header' objects, not HDU objects!
    #
    if fname == None:
        print 'No data files for creating FITS output.'
        raise Exception

    ftemplate = fileutil.openImage(fname,mode='readonly')
    prihdr = pyfits.Header(cards=ftemplate['PRIMARY'].header.ascard.copy())
    del prihdr['pcount']
    del prihdr['gcount']

    if fname.find('.fits') > 0 and len(ftemplate) > 1:

        # Setup which keyword we will use to select each
        # extension...
        _extkey = 'EXTNAME'

        #
        # Now, extract the headers necessary for output (as copies)
        # 1. Find the SCI extension in the template image
        # 2. Make a COPY of the extension header for use in new output file
        extnum = fileutil.findKeywordExtn(ftemplate,_extkey,extlist[0])
        scihdr = pyfits.Header(cards=ftemplate[extnum].header.ascard.copy())

        extnum = fileutil.findKeywordExtn(ftemplate,_extkey,extlist[1])
        errhdr = pyfits.Header(cards=ftemplate[extnum].header.ascard.copy())

        extnum = fileutil.findKeywordExtn(ftemplate,_extkey,extlist[2])
        dqhdr = pyfits.Header(cards=ftemplate[extnum].header.ascard.copy())

    else:
        # Create default headers from scratch
        scihdr = None
        errhdr = None
        dqhdr = None

    ftemplate.close()
    del ftemplate

    # Now, safeguard against having BSCALE and BZERO
    try:
        del scihdr['bscale']
        del scihdr['bzero']
        del errhdr['bscale']
        del errhdr['bzero']
        del dqhdr['bscale']
        del dqhdr['bzero']
    except:
        # If these don't work, they didn't exist to start with...
        pass


    # At this point, check errhdr and dqhdr to make sure they
    # have all the requisite keywords (as listed in updateDTHKeywords).
    # Simply copy them from scihdr if they don't exist...
    if errhdr != None and dqhdr != None:
        for keyword in DTH_KEYWORDS:
            if not errhdr.has_key(keyword):
                errhdr.update(keyword,scihdr[keyword])
            if not dqhdr.has_key(keyword):
                dqhdr.update(keyword,scihdr[keyword])

    return prihdr,scihdr,errhdr,dqhdr
