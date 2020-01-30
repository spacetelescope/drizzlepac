#! /bin/csh

anaconda3

set catcorr=hst_10265_01_acs_wfc_f606w_daophot_corr.cat
set regcorr=${catcorr:r}.reg
cat >! ${regcorr} <<EOF
# Region file format: DS9 version 4.1
# Filename: ${regcorr}
global color=green width=1 font="helvetica 10" select=1 highlite=1 dash=0 fixed=0 edit=0 move=0 delete=0 include=1 source=1
fk5
EOF
tail -n +3 ${catcorr} | awk '{printf "ellipse(%s,%s,0.2\",0.2\",0.0) # tag={%s}\n", $4,$5,$1;}' >> ${regcorr}

set catname=hst_10265_01_acs_wfc_f606w_daophot.cat
set regname=${catname:r}.reg
cat >! ${regname} <<EOF
# Region file format: DS9 version 4.1
# Filename: ${regname}
global color=red width=1 font="helvetica 10" select=1 highlite=1 dash=0 fixed=0 edit=0 move=0 delete=0 include=1 source=1
fk5
EOF
tail -n +3 ${catname} | awk '{printf "ellipse(%s,%s,0.2\",0.2\",0.0) # tag={%s}\n", $4,$5,$1;}' >> ${regname}

set scale=98.0
ds9 \
	hst_10265_01_acs_wfc_f606w_corr.fits \
	-scale mode ${scale} \
	-regions ${regname} \
	hst_10265_01_acs_wfc_f606w_corr.fits \
	-scale mode ${scale} \
	-regions ${regcorr} \
	&

ds9 \
	hst_10265_01_acs_wfc_f606w_drz.fits \
	-scale mode ${scale} \
	-regions ${regname} \
	hst_10265_01_acs_wfc_f606w_drz.fits \
	-scale mode ${scale} \
	-regions ${regcorr} \
	&
