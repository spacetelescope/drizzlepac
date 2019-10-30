"""
This module manages the creation of the output image FITS file.

:Authors: Warren Hack

:License: :doc:`LICENSE`

"""
from astropy.io import fits
from stsci.tools import fileutil, readgeis, logutil

from . import wcs_functions
from . import version
from . import updatehdr

from fitsblender import blendheaders

yes = True
no = False

RESERVED_KEYS = ['NAXIS', 'BITPIX', 'DATE', 'IRAF-TLM', 'XTENSION', 'EXTNAME',' EXTVER']

EXTLIST = ('SCI', 'WHT', 'CTX')

WCS_KEYWORDS = ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'CRPIX1',
'CRPIX2', 'CRVAL1', 'CRVAL2', 'CTYPE1', 'CTYPE2', 'WCSNAME']

# fits.CompImageHDU() crashes with default arguments.
# Instead check that fits module has *attribute* 'CompImageHDU':
PYFITS_COMPRESSION = hasattr(fits, 'CompImageHDU')

# Set up dictionary of default keywords to be written out to the header
# of the output drizzle image using writeDrizKeywords()
DRIZ_KEYWORDS = {
                'VER':{'value':"",'comment':'Drizzle, task version'},
                'GEOM':{'value':"wcs",'comment':'Drizzle, source of geometric information'},
                'DATA':{'value':"",'comment':'Drizzle, input data image'},
                'DEXP':{'value':"",'comment':'Drizzle, input image exposure time (s)'},
                'OUDA':{'value':"",'comment':'Drizzle, output data image'},
                'OUWE':{'value':"",'comment':'Drizzle, output weighting image'},
                'OUCO':{'value':"",'comment':'Drizzle, output context image'},
                'MASK':{'value':"",'comment':'Drizzle, input weighting image'},
                'WTSC':{'value':"",'comment':'Drizzle, weighting factor for input image'},
                'KERN':{'value':"",'comment':'Drizzle, form of weight distribution kernel'},
                'PIXF':{'value':"1.0",'comment':'Drizzle, linear size of drop'},
                'COEF':{'value':"SIP",'comment':'Drizzle, source of coefficients'},
                'OUUN':{'value':"cps",'comment':'Drizzle, units of output image - counts or cps'},
                'FVAL':{'value':"INDEF",'comment':'Drizzle, fill value for zero weight output pix'},
                'WKEY':{'value':"",'comment':'Input image WCS Version used'}
                }

log = logutil.create_logger(__name__, level=logutil.logging.NOTSET)


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
        The object 'plist' must contain at least the following members::

            plist['output']   - name of output FITS image (for SCI)
            plist['outnx']    - size of X axis for output array
            plist['outny']    - size of Y axis for output array

        If 'single=yes', then 'plist' also needs to contain::

            plist['outsingle']
            plist['outsweight']
            plist['outscontext']

        If 'blot=yes', then 'plist' also needs::

            plist['data']
            plist['blotImage']
            plist['blotnx'],plist['blotny']

        If 'build' is set to 'no', then each extension/array must be
        in separate FITS objects.  This would also require::

          plist['outdata']    - name of output SCI FITS image
          plist['outweight']  - name of output WHT FITS image
          plist['outcontext'] - name of output CTX FITS image

        Optionally, the overall exposure time information can be passed as::

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
        self.blot = blot

        if PYFITS_COMPRESSION and 'compress' in input_pars:
            self.compress = input_pars['compress'] # Control creation of compressed FITS files
        else:
            self.compress = False

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
            self.texptime = plist[0]['exptime']
            self.expstart = plist[0]['expstart']
            self.expend = plist[0]['expend']
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
            self.outblot_key = 'blotImage'

        if not self.build or single:
            self.output = _outdata

        self.outdata = _outdata
        self.outweight = _outweight
        self.outcontext = _outcontext

    def set_bunit(self,bunit):
        """
        Method used to update the value of the bunit attribute.
        """
        self.bunit = bunit

    def set_units(self,units):
        """
        Method used to record what units were specified by the user for the output product.
        """
        self.units = units

    def writeFITS(self, template, sciarr, whtarr, ctxarr=None,
                versions=None, overwrite=yes, blend=True, virtual=False):
        """
        Generate PyFITS objects for each output extension
        using the file given by 'template' for populating
        headers.

        The arrays will have the size specified by 'shape'.
        """
        if not isinstance(template, list):
            template = [template]

        if fileutil.findFile(self.output):
            if overwrite:
                log.info('Deleting previous output product: %s' % self.output)
                fileutil.removeFile(self.output)

            else:
                log.warning('Output file %s already exists and overwrite not '
                            'specified!' % self.output)
                log.error('Quitting... Please remove before resuming '
                          'operations.')
                raise IOError

        # initialize output value for this method
        outputFITS = {}
        # Default value for NEXTEND when 'build'== True
        nextend = 3
        if not self.build:
            nextend = 0
            if self.outweight:
                if overwrite:
                    if fileutil.findFile(self.outweight):
                        log.info('Deleting previous output WHT product: %s' %
                                 self.outweight)
                    fileutil.removeFile(self.outweight)
                else:
                    log.warning('Output file %s already exists and overwrite '
                                'not specified!' % self.outweight)
                    log.error('Quitting... Please remove before resuming '
                              'operations.')
                    raise IOError


            if self.outcontext:
                if overwrite:
                    if fileutil.findFile(self.outcontext):
                        log.info('Deleting previous output CTX product: %s' %
                                 self.outcontext)
                    fileutil.removeFile(self.outcontext)
                else:
                    log.warning('Output file %s already exists and overwrite '
                                'not specified!' % self.outcontext)
                    log.error('Quitting... Please remove before resuming '
                              'operations.')
                    raise IOError


        # Get default headers from multi-extension FITS file
        # If only writing out single drizzle product, blending needs to be
        # forced off as there is only 1 input to report, no blending needed
        if self.single:
            blend=False

        # If input data is not in MEF FITS format, it will return 'None'
        # and those headers will have to be generated from drizzle output
        # file FITS headers.
        # NOTE: These are HEADER objects, not HDUs
        #prihdr,scihdr,errhdr,dqhdr = getTemplates(template)
        self.fullhdrs, intab = getTemplates(template, blend=False)

        newhdrs, newtab = getTemplates(template,blend=blend)
        if newtab is not None: nextend += 1 # account for new table extn

        prihdr = newhdrs[0]
        scihdr = newhdrs[1]
        errhdr = newhdrs[2]
        dqhdr = newhdrs[3]

        # Setup primary header as an HDU ready for appending to output FITS file
        prihdu = fits.PrimaryHDU(header=prihdr, data=None)

        # Start by updating PRIMARY header keywords...
        prihdu.header.set('EXTEND', value=True, after='NAXIS')
        prihdu.header['NEXTEND'] = nextend
        prihdu.header['FILENAME'] = self.output
        prihdu.header['PROD_VER'] = 'DrizzlePac {}'.format(version.__version__)

        # Update the ROOTNAME with the new value as well
        _indx = self.output.find('_drz')
        if _indx < 0:
            rootname_val = self.output
        else:
            rootname_val = self.output[:_indx]
        prihdu.header['ROOTNAME'] = rootname_val


        # Get the total exposure time for the image
        # If not calculated by PyDrizzle and passed through
        # the pardict, then leave value from the template image.
        if self.texptime:
            prihdu.header['EXPTIME'] = self.texptime
            prihdu.header.set('TEXPTIME', value=self.texptime, after='EXPTIME')
            prihdu.header['EXPSTART'] = self.expstart
            prihdu.header['EXPEND'] = self.expend

        #Update ASN_MTYPE to reflect the fact that this is a product
        # Currently hard-wired to always output 'PROD-DTH' as MTYPE
        prihdu.header['ASN_MTYP'] = 'PROD-DTH'

        # Update DITHCORR calibration keyword if present
        # Remove when we can modify FITS headers in place...
        if 'DRIZCORR' in prihdu.header:
            prihdu.header['DRIZCORR'] = 'COMPLETE'
        if 'DITHCORR' in prihdu.header:
            prihdu.header['DITHCORR'] = 'COMPLETE'

        prihdu.header['NDRIZIM'] =(len(self.parlist),
                                   'Drizzle, No. images drizzled onto output')

        # Only a subset of these keywords makes sense for the new WCS based
        # transformations. They need to be reviewed to decide what to keep
        # and what to leave out.
        if not self.blot:
            self.addDrizKeywords(prihdu.header,versions)

        if scihdr:
            try:
                del scihdr['OBJECT']
            except KeyError:
                pass

            if 'CCDCHIP' in scihdr: scihdr['CCDCHIP'] = '-999'
            if 'NCOMBINE' in scihdr:
                scihdr['NCOMBINE'] = self.parlist[0]['nimages']

            # If BUNIT keyword was found and reset, then
            bunit_last_kw = self.find_kwupdate_location(scihdr,'bunit')
            if self.bunit is not None:
                comment_str = "Units of science product"
                if self.bunit.lower()[:5] == 'count':
                    comment_str = "counts * gain = electrons"
                scihdr.set('BUNIT', value=self.bunit,
                           comment=comment_str,
                           after=bunit_last_kw)
            else:
                # check to see whether to update already present BUNIT comment
                if 'bunit' in scihdr and scihdr['bunit'].lower()[:5] == 'count':
                    comment_str = "counts * gain = electrons"
                    scihdr.set('BUNIT', value=scihdr['bunit'],
                               comment=comment_str,
                               after=bunit_last_kw)

            # Add WCS keywords to SCI header
            if self.wcs:
                pre_wcs_kw = self.find_kwupdate_location(scihdr,'CD1_1')
                addWCSKeywords(self.wcs,scihdr,blot=self.blot,
                                single=self.single, after=pre_wcs_kw)
                # Recompute this after removing distortion kws
                pre_wcs_kw = self.find_kwupdate_location(scihdr,'CD1_1')

        ##########
        # Now, build the output file
        ##########
        if self.build:
            print('-Generating multi-extension output file: ',self.output)
            fo = fits.HDUList()

            # Add primary header to output file...
            fo.append(prihdu)

            if self.single and self.compress:
                hdu = fits.CompImageHDU(data=sciarr, header=scihdr, name=EXTLIST[0])
            else:
                hdu = fits.ImageHDU(data=sciarr, header=scihdr, name=EXTLIST[0])
            last_kw = self.find_kwupdate_location(scihdr,'EXTNAME')
            hdu.header.set('EXTNAME', value='SCI', after=last_kw)
            hdu.header.set('EXTVER', value=1, after='EXTNAME')
            fo.append(hdu)

            # Build WHT extension here, if requested...
            if errhdr:
                errhdr['CCDCHIP'] = '-999'

            if self.single and self.compress:
                hdu = fits.CompImageHDU(data=whtarr, header=errhdr, name=EXTLIST[1])
            else:
                hdu = fits.ImageHDU(data=whtarr, header=errhdr, name=EXTLIST[1])
            last_kw = self.find_kwupdate_location(errhdr,'EXTNAME')
            hdu.header.set('EXTNAME', value='WHT', after=last_kw)
            hdu.header.set('EXTVER', value=1, after='EXTNAME')
            if self.wcs:
                pre_wcs_kw = self.find_kwupdate_location(hdu.header,'CD1_1')
                # Update WCS Keywords based on PyDrizzle product's value
                # since 'drizzle' itself doesn't update that keyword.
                addWCSKeywords(self.wcs,hdu.header,blot=self.blot,
                               single=self.single, after=pre_wcs_kw)
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

            if self.single and self.compress:
                hdu = fits.CompImageHDU(data=_ctxarr, header=dqhdr, name=EXTLIST[2])
            else:
                hdu = fits.ImageHDU(data=_ctxarr, header=dqhdr, name=EXTLIST[2])
            last_kw = self.find_kwupdate_location(dqhdr,'EXTNAME')
            hdu.header.set('EXTNAME', value='CTX', after=last_kw)
            hdu.header.set('EXTVER', value=1, after='EXTNAME')

            if self.wcs:
                pre_wcs_kw = self.find_kwupdate_location(hdu.header,'CD1_1')
                # Update WCS Keywords based on PyDrizzle product's value
                # since 'drizzle' itself doesn't update that keyword.
                addWCSKeywords(self.wcs,hdu.header,blot=self.blot,
                               single=self.single, after=pre_wcs_kw)
            fo.append(hdu)

            # remove all alternate WCS solutions from headers of this product
            wcs_functions.removeAllAltWCS(fo,[1])

            # add table of combined header keyword values to FITS file
            if newtab is not None:
                fo.append(newtab)

            if not virtual:
                print('Writing out to disk:',self.output)
                # write out file to disk
                fo.writeto(self.output)
                fo.close()
                del fo, hdu
                fo = None
            # End 'if not virtual'
            outputFITS[self.output]= fo

        else:
            print('-Generating simple FITS output: %s' % self.outdata)

            fo = fits.HDUList()
            hdu_header = prihdu.header.copy()
            del hdu_header['nextend']

            # Append remaining unique header keywords from template DQ
            # header to Primary header...
            if scihdr:
                for _card in scihdr.cards:
                    if _card.keyword not in RESERVED_KEYS and _card.keyword not in hdu_header:
                        hdu_header.append(_card)
            for kw in ['PCOUNT', 'GCOUNT']:
                try:
                    del kw
                except KeyError:
                    pass
            hdu_header['filename'] = self.outdata

            if self.compress:
                hdu = fits.CompImageHDU(data=sciarr, header=hdu_header)
                wcs_ext = [1]
            else:
                hdu = fits.ImageHDU(data=sciarr, header=hdu_header)
                wcs_ext = [0]

            # explicitly set EXTEND to FALSE for simple FITS files.
            dim = len(sciarr.shape)
            hdu.header.set('extend', value=False, after='NAXIS%s'%dim)

            # Add primary header to output file...
            fo.append(hdu)

            # remove all alternate WCS solutions from headers of this product
            logutil.logging.disable(logutil.logging.INFO)
            wcs_functions.removeAllAltWCS(fo,wcs_ext)
            logutil.logging.disable(logutil.logging.NOTSET)

            # add table of combined header keyword values to FITS file
            if newtab is not None:
                fo.append(newtab)

            if not virtual or "single_sci" in self.outdata:
                print('Writing out image to disk:',self.outdata)
                # write out file to disk
                fo.writeto(self.outdata)
                del hdu
                if "single_sci" not in self.outdata:
                    del fo
                    fo = None
            # End 'if not virtual'
            outputFITS[self.outdata]= fo

            if self.outweight and whtarr is not None:
                # We need to build new PyFITS objects for each WHT array
                fwht = fits.HDUList()

                if errhdr:
                    errhdr['CCDCHIP'] = '-999'

                if self.compress:
                    hdu = fits.CompImageHDU(data=whtarr, header=prihdu.header)
                else:
                    hdu = fits.ImageHDU(data=whtarr, header=prihdu.header)
                # Append remaining unique header keywords from template DQ
                # header to Primary header...
                if errhdr:
                    for _card in errhdr.cards:
                        if _card.keyword not in RESERVED_KEYS and _card.keyword not in hdu.header:
                            hdu.header.append(_card)
                hdu.header['filename'] = self.outweight
                hdu.header['CCDCHIP'] = '-999'
                if self.wcs:
                    pre_wcs_kw = self.find_kwupdate_location(hdu.header,'CD1_1')
                    # Update WCS Keywords based on PyDrizzle product's value
                    # since 'drizzle' itself doesn't update that keyword.
                    addWCSKeywords(self.wcs,hdu.header, blot=self.blot,
                                   single=self.single, after=pre_wcs_kw)

                # Add primary header to output file...
                fwht.append(hdu)
                # remove all alternate WCS solutions from headers of this product
                wcs_functions.removeAllAltWCS(fwht,wcs_ext)

                if not virtual:
                    print('Writing out image to disk:',self.outweight)
                    fwht.writeto(self.outweight)
                    del fwht,hdu
                    fwht = None
                # End 'if not virtual'
                outputFITS[self.outweight]= fwht

            # If a context image was specified, build a PyFITS object
            # for it as well...
            if self.outcontext and ctxarr is not None:
                fctx = fits.HDUList()

                # If there is only 1 plane, write it out as a 2-D extension
                if ctxarr.shape[0] == 1:
                    _ctxarr = ctxarr[0]
                else:
                    _ctxarr = ctxarr

                if self.compress:
                    hdu = fits.CompImageHDU(data=_ctxarr, header=prihdu.header)
                else:
                    hdu = fits.ImageHDU(data=_ctxarr, header=prihdu.header)
                # Append remaining unique header keywords from template DQ
                # header to Primary header...
                if dqhdr:
                    for _card in dqhdr.cards:
                        if ( (_card.keyword not in RESERVED_KEYS) and
                             _card.keyword not in hdu.header):
                            hdu.header.append(_card)
                hdu.header['filename'] = self.outcontext
                if self.wcs:
                    pre_wcs_kw = self.find_kwupdate_location(hdu.header,'CD1_1')
                    # Update WCS Keywords based on PyDrizzle product's value
                    # since 'drizzle' itself doesn't update that keyword.
                    addWCSKeywords(self.wcs,hdu.header, blot=self.blot,
                                   single=self.single, after=pre_wcs_kw)

                fctx.append(hdu)
                # remove all alternate WCS solutions from headers of this product
                wcs_functions.removeAllAltWCS(fctx,wcs_ext)
                if not virtual:
                    print('Writing out image to disk:',self.outcontext)
                    fctx.writeto(self.outcontext)
                    del fctx,hdu
                    fctx = None
                # End 'if not virtual'

                outputFITS[self.outcontext]= fctx

        return outputFITS

    def find_kwupdate_location(self,hdr,keyword):
        """
        Find the last keyword in the output header that comes before the new
        keyword in the original, full input headers.
        This will rely on the original ordering of keywords from the original input
        files in order to place the updated keyword in the correct location in case
        the keyword was removed from the output header prior to calling this method.
        """
        # start by looping through the full templates
        kw_list = None
        last_kw = None
        for extn in self.fullhdrs:
            if keyword in extn:
                #indx = extn.ascard.index_of(keyword)
                indx = extn.index(keyword)
                kw_list = list(extn.keys())[:indx]
                break
        if kw_list:
            # find which keyword from this list exists in header to be updated
            for kw in kw_list[::-1]:
                if kw in hdr:
                    last_kw = kw
                    break
        # determine new value for the last keyword found before the HISTORY kws
        if last_kw is None:
            hdrkeys = list(hdr.keys())
            i = -1
            last_kw = hdrkeys[i]
            while last_kw == 'HISTORY':
                i -= 1
                last_kw = hdrkeys[i]

        return last_kw


    def addDrizKeywords(self,hdr,versions):
        """ Add drizzle parameter keywords to header. """

        # Extract some global information for the keywords
        _geom = 'User parameters'

        _imgnum = 0
        for pl in self.parlist:

            # Start by building up the keyword prefix based
            # on the image number for the chip
            #_keyprefix = 'D%03d'%_imgnum
            _imgnum += 1

            drizdict = DRIZ_KEYWORDS.copy()
            # Update drizdict with current values
            drizdict['VER']['value'] = pl['driz_version'][:44]
            drizdict['DATA']['value'] = pl['data'][:64]
            drizdict['DEXP']['value'] = pl['exptime']
            drizdict['OUDA']['value'] = pl['outFinal'][:64]
            drizdict['OUWE']['value'] = pl['outWeight'][:64]
            if pl['outContext'] is None:
                outcontext = ""
            else:
                outcontext = pl['outContext'][:64]
            drizdict['OUCO']['value'] = outcontext
            if self.single:
                drizdict['MASK']['value'] = pl['singleDrizMask'][:64]
            else:
                drizdict['MASK']['value'] = pl['finalMask'][:64]

            # Process the values of WT_SCL to be consistent with
            # what IRAF Drizzle would output
            if 'wt_scl_val' in pl:
                _wtscl = pl['wt_scl_val']
            else:
                if pl['wt_scl'] == 'exptime': _wtscl = pl['exptime']
                elif pl['wt_scl'] == 'expsq': _wtscl = pl['exptime']*pl['exptime']
                else: _wtscl = pl['wt_scl']

            drizdict['WTSC']['value'] = _wtscl
            drizdict['KERN']['value'] = pl['kernel']
            drizdict['PIXF']['value'] = pl['pixfrac']
            drizdict['OUUN']['value'] = self.units
            if pl['fillval'] is None:
                _fillval = 'INDEF'
            else:
                _fillval = pl['fillval']
            drizdict['FVAL']['value'] = _fillval
            drizdict['WKEY']['value'] = pl['driz_wcskey']

            drizdict['SCAL'] = {'value':pl['scale'],'comment':'Drizzle, pixel size (arcsec) of output image'}
            drizdict['ISCL'] = {'value':pl['idcscale'],'comment':'Drizzle, default IDCTAB pixel size(arcsec)'}

            # Now update header with values
            writeDrizKeywords(hdr,_imgnum,drizdict)
            del drizdict

        # Add version information as HISTORY cards to the header
        if versions is not None:
            ver_str = "AstroDrizzle processing performed using: "
            hdr.add_history(ver_str)
            for k in versions.keys():
                ver_str = '    '+str(k)+' Version '+str(versions[k])
                hdr.add_history(ver_str)


def cleanTemplates(scihdr,errhdr,dqhdr):

    # Now, safeguard against having BSCALE and BZERO
    for kw in ['BSCALE', 'BZERO']:
        try:
            del scihdr[kw]
        except KeyError:
            pass
        try:
            del errhdr[kw]
        except KeyError:
            pass
        try:
            del dqhdr[kw]
        except KeyError:
            pass

    # At this point, check errhdr and dqhdr to make sure they
    # have all the requisite keywords (as listed in updateDTHKeywords).
    # Simply copy them from scihdr if they don't exist...
    if errhdr is not None and dqhdr is not None:
        for keyword in WCS_KEYWORDS:
            if keyword in scihdr:
                if keyword not in errhdr:
                    errhdr[keyword] = scihdr[keyword]
                if keyword not in dqhdr:
                    dqhdr[keyword]= scihdr[keyword]

def getTemplates(fnames, blend=True):
    """ Process all headers to produce a set of combined headers
        that follows the rules defined by each instrument.

    """
    if not blend:
        newhdrs =  blendheaders.getSingleTemplate(fnames[0])
        newtab = None
    else:
        # apply rules to create final version of headers, plus table
        newhdrs, newtab = blendheaders.get_blended_headers(inputs=fnames)

    cleanTemplates(newhdrs[1],newhdrs[2],newhdrs[3])

    return newhdrs, newtab

def addWCSKeywords(wcs,hdr,blot=False,single=False,after=None):
    """ Update input header 'hdr' with WCS keywords.
    """
    wname = wcs.wcs.name
    wtype = updatehdr.interpret_wcsname_type(wname)

    # Update WCS Keywords based on PyDrizzle product's value
    # since 'drizzle' itself doesn't update that keyword.
    hdr['WCSNAME'] = wname
    hdr['WCSTYPE'] = wtype
    hdr.set('VAFACTOR', value=1.0, after=after)
    hdr.set('ORIENTAT', value=wcs.orientat, after=after)

    # Use of 'after' not needed if these keywords already exist in the header
    if after in WCS_KEYWORDS:
        after = None

    if 'CTYPE1' not in hdr:
        hdr.set('CTYPE2', value=wcs.wcs.ctype[1], after=after)
        hdr.set('CTYPE1', value=wcs.wcs.ctype[0], after=after)
    hdr.set('CRPIX2', value=wcs.wcs.crpix[1], after=after)
    hdr.set('CRPIX1', value=wcs.wcs.crpix[0], after=after)
    hdr.set('CRVAL2', value=wcs.wcs.crval[1], after=after)
    hdr.set('CRVAL1', value=wcs.wcs.crval[0], after=after)
    hdr.set('CD2_2', value=wcs.wcs.cd[1][1], after=after)
    hdr.set('CD2_1', value=wcs.wcs.cd[1][0], after=after)
    hdr.set('CD1_2', value=wcs.wcs.cd[0][1], after=after)
    hdr.set('CD1_1', value=wcs.wcs.cd[0][0], after=after)

    # delete distortion model related keywords
    deleteDistortionKeywords(hdr)

    if not blot:
        blendheaders.remove_distortion_keywords(hdr)

def deleteDistortionKeywords(hdr):
    """ Delete distortion related keywords from output drizzle science header
        since the drizzled image should have no remaining distortion.
    """
    dist_kws = ['D2IMERR1','D2IMERR2','D2IMDIS1','D2IMDIS2','D2IM1.*','D2IM2.*','D2IMEXT']
    for kw in dist_kws:
        try:
            del hdr[kw]
        except KeyError as e:
            pass


def writeSingleFITS(data,wcs,output,template,clobber=True,verbose=True):
    """ Write out a simple FITS file given a numpy array and the name of another
    FITS file to use as a template for the output image header.
    """
    outname,outextn = fileutil.parseFilename(output)
    outextname,outextver = fileutil.parseExtn(outextn)

    if fileutil.findFile(outname):
        if clobber:
            log.info('Deleting previous output product: %s' % outname)
            fileutil.removeFile(outname)

        else:
            log.warning('Output file %s already exists and overwrite not '
                        'specified!' % outname)
            log.error('Quitting... Please remove before resuming operations.')
            raise IOError

    # Now update WCS keywords with values from provided WCS
    if hasattr(wcs.sip,'a_order'):
        siphdr = True
    else:
        siphdr = False
    wcshdr = wcs.wcs2header(sip2hdr=siphdr)

    if template is not None:
        # Get default headers from multi-extension FITS file
        # If input data is not in MEF FITS format, it will return 'None'
        # NOTE: These are HEADER objects, not HDUs
        (prihdr,scihdr,errhdr,dqhdr),newtab = getTemplates(template,EXTLIST)

        if scihdr is None:
            scihdr = fits.Header()
            indx = 0
            for c in prihdr.cards:
                if c.keyword not in ['INHERIT','EXPNAME']: indx += 1
                else: break
            for i in range(indx,len(prihdr)):
                scihdr.append(prihdr.cards[i])
            for i in range(indx, len(prihdr)):
                del prihdr[indx]
    else:
        scihdr = fits.Header()
        prihdr = fits.Header()
        # Start by updating PRIMARY header keywords...
        prihdr.set('EXTEND', value=True, after='NAXIS')
        prihdr['FILENAME'] = outname

    if outextname == '':
        outextname = 'sci'
    if outextver == 0: outextver = 1
    scihdr['EXTNAME'] = outextname.upper()
    scihdr['EXTVER'] = outextver

    for card in wcshdr.cards:
        scihdr[card.keyword] = (card.value, card.comment)

    # Create PyFITS HDUList for all extensions
    outhdu = fits.HDUList()
    # Setup primary header as an HDU ready for appending to output FITS file
    prihdu = fits.PrimaryHDU(header=prihdr)
    scihdu = fits.ImageHDU(header=scihdr,data=data)

    outhdu.append(prihdu)
    outhdu.append(scihdu)
    outhdu.writeto(outname)

    if verbose:
        print('Created output image: %s' % outname)

def writeDrizKeywords(hdr,imgnum,drizdict):
    """ Write basic drizzle-related keywords out to image header as a record
        of the processing performed to create the image

        The dictionary 'drizdict' will contain the keywords and values to be
        written out to the header.
    """
    _keyprefix = 'D%03d'%imgnum

    for key in drizdict:
        val = drizdict[key]['value']
        if val is None: val = ""
        comment = drizdict[key]['comment']
        if comment is None: comment = ""
        hdr[_keyprefix+key] = (val, drizdict[key]['comment'])
