include <iraf77.h>

# TIMOTP -- Temporary replacement for uimotp

procedure timotp (f77nam, listp, istat)

%	character*(*)	f77nam
pointer	listp		# o: image template pointer
int	istat		# o: return status
#--
pointer	sp, sppnam
pointer	imtopen()

begin
	# Convert fortran to spp string

	call smark (sp)
	call salloc (sppnam, SZ_PATHNAME, TY_CHAR)

	call f77upk (f77nam, Memc[sppnam], SZ_PATHNAME)

	# Call spp routine with error checking

	istat = ER_OK
	iferr (listp = imtopen (Memc[sppnam])) {
	    istat = ER_IMTOPN
	    listp = NULL
	}

	call sfree (sp)
end

# TIMXTP -- Temporary replacement for uimxtp

procedure timxtp (listp, f77nam, istat)

pointer	listp		# i: image template pointer
%	character*(*)	f77nam
int	istat		# o: return status
#--
pointer	sp, sppnam
int	imtgetim()

begin
	# Allocate dynamic memory for temporary string

	call smark (sp)
	call salloc (sppnam, SZ_PATHNAME, TY_CHAR)

	# Call spp routine with error checking

	iferr (istat = imtgetim (listp, Memc[sppnam], SZ_PATHNAME)) {
	    istat = ER_IMTGETNEXT
	    Memc[sppnam] = EOS
	} else if (istat == EOF) {
	    istat = ER_EOF
	} else {
	    istat = ER_OK
	}

	# Convert spp to fortran string

	call f77pak (Memc[sppnam], f77nam, SZ_PATHNAME)
	call sfree (sp)
end

# TIMCTP -- Temporary replacement for oimtcp

procedure timctp (listp, istat)

pointer	listp		# i: image template pointer
int	istat		# o: return status
#--

begin
	# Call spp routine with error checking

	istat = ER_OK
	iferr (call imtclose (listp))
	    istat = ER_IMTCLST

end
