;+
; NAME:
;  megsb_image_to_spectrum
;
; PURPOSE:
;  Create a constant wavelength scale spectrum from a MEGS-B image
;
; CATEGORY:
;  Level 1 MEGS-B
;
; CALLING SEQUENCE:
;  megsb_image_to_spectrum, img, spwave, spflux,
;      mask=mask, err=err, good=good, imgwave=imgwave, sperr=sperr,
;     linefile='thelinefile'
;
; INPUTS:
;  img: a 2048x1024 MEGS-B image (ff corrected)
;
; OPTIONAL INPUTS:
;  mask: if indefined, a default mask is used and applied to the img
;  err: a relative uncertainty in each pixel of img
;  good: a bytarray of 0 (bad) or 1 (good) pixels in img
;  imgwave: a 2048x1024 wavelength array for img (no fitting)
;
; KEYWORD PARAMETERS:
;  linefile: a file name to use for fitting lines
;   (if you don't know, don't use it)
;
; OUTPUTS:
;  spwave: the uniform wavelength array
;  spflux: the flux array for each wavelength bin in spwave
;
; OPTIONAL OUTPUTS:
;  sperr: spectrum flux uncertainty propagated from err
;  mask & imgwave will contain valid values if they point to undefined variables
;
; COMMON BLOCKS:
;  megsb_image_to_spectrum_cal:
;   default_mask - the default mask from a saveset
;   default_imgwave - a default wavelength scale that should work for
;                     center pointing
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;  megsb_image_to_spectrum, image, wave, flux
;   ; This example uses the default wavelength scale.
;
; MODIFICATION HISTORY:
;  10/06/06 DLW Original file creation
;
;-

PRO megsb_image_to_spectrum, imgfull, spwave, spflux, $
	imgmask=imgmask, imgacc=imgacc, imgprec=imgprec, imgwgd=imgwgd, $
	imgwave=imgwave, spextra1=spextra1, spextra2=spextra2, linefile=linefile, $
	spaccout=spaccout, spprecout=spprecout, $
	spextra1out=spextra1out, spextra2out=spextra2out

	; imgfull is a 2048x1024 ff corrected image from megs-b
	; imgmask is a 2048x1024 mask to apply to the image (if undefined, a default mask is used)
	; imgacc is a 2048x1024 float array of relative uncertainties (accuracy) to propagate
	; imgprec is a 2048x1024 float array of relative uncertainties (precision) to propagate
	; imgwgd is a 2048x1024 bytarr of 0 (bad) or 1 (good) indicating where imgfull is good
	; imgwave is an optional 2048x1024 wavelength image to use
	;  (actual wavelengths are returned if not specified)
	; spwave is the one-dimensional wavelength array used in regridding
	; spflux is the output spectral flux
	; spextra1 is an optional array to regrid (regrid is spextra1out)
	; spextra2 is an optional array to regrid (regrid is spextra2out)
	; linefile=linefile allows the user to specify a file of lines to use for fitting (overriding any imgwave specified)
	; spaccout is spectrum accuracy propagated
	; spprecout is spectrum precision propagated

	COMMON megsb_image_to_spectrum_cal, default_mask, default_imgwave

;	@megs_defines

        @megs_wave_defines
        
	IF SIZE(default_mask,/type) EQ 0 THEN BEGIN
		;setup common block
		   ;restore,'megsb_default_mask.sav'
           default_mask = imgmask ; 10/2/19 define default mask for rocket
           default_imgwave = imgwave
;		RESTORE, GETENV('eve_cal_data') + $
;			'/megsb_default_mask.sav' ;restores default_mask variable
;		IF is_rocket() EQ 0 THEN BEGIN
;			; flight
;			RESTORE, GETENV('eve_cal_data') + $
;				'/megsb_wave.sav'     ;restores waveimg ;megsb_wave
;			default_imgwave = TEMPORARY(waveimg) ;assign to default_imgwave variable
;		ENDIF ELSE BEGIN
			; rocket
;			RESTORE,GETENV('eve_cal_data')+'/rkt_36240_megsb_full_wave.sav' ;waveimg
;			default_imgwave = TEMPORARY(waveimg) ;assign to default_imgwave variable
;		ENDELSE

		; provide fill values for missing fit data
		;  fill vertically
		;tvscl,congrid(default_imgwave,1024,512)
		FOR i=0L,2047 DO BEGIN
			x=WHERE(default_imgwave[i,*] LT 1.,n_x,comp=comp,ncomp=ncomp)
			IF n_x GT 0 THEN BEGIN
				default_imgwave[i,0:comp[0]]=default_imgwave[i,comp[0]]
				default_imgwave[i,comp[ncomp-1]:1023]=default_imgwave[i,comp[ncomp-1]]
			ENDIF
		ENDFOR
		;tvscl,congrid(default_imgwave,1024,512)

	ENDIF

	; is a mask provided?
	IF SIZE(imgmask,/type) EQ 0 THEN imgmask = default_mask

	; COMMENTED OUT 05/21/08 DLW
	;megsb_tuned_mask, imgfull, imgmask, tunedmask
	; imgmask has all the warts
	; tunedmask is a trimmed-down imgmask
	;imgmask = imgmask AND tunedmask

	; get the stripe image
	stripeimg = make_megsb_stripe( imgfull, imgmask )

	; try to force rocket environment to use default wavelengths only
;	IF is_rocket() THEN BEGIN
;		linefile=''
;		tmp=linefile ;undefine it
;		imgwave = default_imgwave
;	ENDIF

	IF SIZE(imgwave,/type) EQ 0 THEN BEGIN
		;imgwave is undefined, may need to do the full fit
		IF SIZE(linefile,/type) EQ 0 THEN BEGIN
			; use default_imgwave
			PRINT,'MEGSB_IMAGE_TO_SPECTRUM using default wavelengths'
			thisimgwave = make_megsb_stripe( default_imgwave, imgmask )
			; ERRORS ARE CURRENTLY IGNORED
			;make_megsb_wave_img, stripeimg, thisimgwave, spwave, spflux ;, err=err
			make_megsb_wave_img, stripeimg, thisimgwave, spwave, spflux, $
				imgmask=imgmask, imgacc=imgacc, imgprec=imgprec, imgwgd=imgwgd, $
				imgwave=imgwave, spextra1=spextra1, spextra2=spextra2, $
				spaccout=spaccout, spprecout=spprecout, $
				spextra1out=spextra1out, spextra2out=spextra2out
		ENDIF ELSE BEGIN
			;do the full fit with a specific linefile (the works!)
			PRINT,'MEGSB_IMAGE_TO_SPECTRUM fitting wavelengths'
			;make_megsb_wave_img, stripeimg, imgwave, spwave, spflux, file=linefile
			make_megsb_wave_img, stripeimg, imgwave, spwave, spflux, file=linefile, $
				imgmask=imgmask, imgacc=imgacc, imgprec=imgprec, imgwgd=imgwgd, $
				imgwave=imgwave, spextra1=spextra1, spextra2=spextra2, $
				spaccout=spaccout, spprecout=spprecout, $
				spextra1out=spextra1out, spextra2out=spextra2out
			;stop
		ENDELSE
	ENDIF ELSE BEGIN
		; get the stripe waveimg if needed
		;print,'MEGSB_IMAGE_TO_SPECTRUM using wavelength argument (not a fit)'
		IF N_ELEMENTS(imgwave) NE N_ELEMENTS(stripeimg) THEN $
			thisimgwave = make_megsb_stripe( imgwave, imgmask ) $
		ELSE thisimgwave = imgwave

		; fill the wavelength array to avoid div by zero
		FOR i=0L,2047 DO BEGIN
			x=WHERE(thisimgwave[i,*] LT 1.,n_x,comp=comp,ncomp=ncomp)
			IF n_x GT 0 AND N_ELEMENTS(thisimgwave[i,*])-1 NE n_x THEN BEGIN
				IF comp[0] GT 0 THEN $
					thisimgwave[i,0:comp[0]]=thisimgwave[i,comp[0]]
				IF ncomp GT 0 AND ncomp-1 LT N_ELEMENTS(thisimgwave[i,*])-1 THEN $
					thisimgwave[i,comp[ncomp-1]:*]=thisimgwave[i,comp[ncomp-1]]
			ENDIF
		ENDFOR
		;***
		;***
		; need to fill thisimgwave horizontally now
		;***
		;***

		; use the input wavelength scale instead of fitting
		n_spec = N_ELEMENTS(stripeimg[0,*])
		n_wave = N_ELEMENTS(stripeimg[*,0])

		IF SIZE(spextra1,/type) EQ 0 THEN spextra1=FLTARR(n_spec,n_wave)>1.
		IF N_ELEMENTS(spextra1) NE N_ELEMENTS(stripeimg) THEN $
			sp_extra1 = make_megsb_stripe( spextra1, imgmask ) $
		ELSE sp_extra1=spextra1
		IF SIZE(spextra2,/type) EQ 0 THEN spextra2=spextra1
		IF N_ELEMENTS(spextra2) NE N_ELEMENTS(stripeimg) THEN $
			sp_extra2 = make_megsb_stripe( spextra2, imgmask ) $
		ELSE sp_extra2=spextra2

		; only get the relevant stripe
		IF SIZE(imgwgd,/type) EQ 0 THEN imgwgd=BYTARR(2048,1024)>1
		IF N_ELEMENTS(imgwgd) NE N_ELEMENTS(stripeimg) THEN $
			thisimgwgd = make_megsb_stripe( imgwgd, imgmask ) $
		ELSE thisimgwgd = BYTARR(n_spec,n_wave)>1b

		; only get the relevant stripe
		IF SIZE(imgacc,/type) EQ 0 THEN imgacc=FLTARR(2048,1024)>1.
		IF N_ELEMENTS(imgacc) NE N_ELEMENTS(stripeimg) THEN $
			thisimgacc = make_megsb_stripe( imgacc, imgmask )
		;else thisimgacc = fltarr(n_spec,n_wave)>1.

		; only get the relevant stripe
		IF SIZE(imgprec,/type) EQ 0 THEN imgprec=thisimgacc
		thisimgprec=imgprec
		IF N_ELEMENTS(imgprec) NE N_ELEMENTS(stripeimg) THEN $
			thisimgprec = make_megsb_stripe( imgprec, imgmask )
		;else thisimgprec = fltarr(n_spec,n_wave)>1.

		spwave = MEGSB_L1_WAVE ;dindgen(3100)*.025d0 + 31.d0
		n_out  = N_ELEMENTS(spwave)
		spflux_out = DBLARR(n_out)

		;help,thisimgprec
		;stop
		;regrid_spectra, thisimgwave>.1, stripeimg > 1.e-14, n_spec, n_wave, thisimgwgd, $
		regrid_spectra, thisimgwave>.1, stripeimg, n_spec, n_wave, thisimgwgd, $
			thisimgacc, thisimgprec, $
			spwave, spflux_out, n_out, wg_out, $
			spprec_out, spacc_out, $
			sp_extra1, sp_extra2, $
			spextra1out, spextra2out, status
		;stop

		; force accuracy to be no less than precision
		spprec_out = spprec_out > 1e-12
		spacc_out = spacc_out > spprec_out
		; create a mask to remove low signal areas that bias the errors high
		IF SIZE(spflux_out,/n_dim) EQ 2 THEN BEGIN
			sp = TOTAL(spflux_out,2) / n_spec
		ENDIF ELSE sp = spflux_out
		im = REBIN(sp>1e-12 ,n_out,n_spec)
		good_flux_mask = (im LT spflux_out) ; limit by the column average

		; precision	### NEED TO REVISIT
;		good_err_mask = (spprec_out LT (spflux_out > 1e-12)) AND (spflux_out GT 0.0) ; byte mask
		good_err_mask = (spflux_out > 1e-12) AND (spflux_out GT 0.0) ; byte mask

		;propagate precision the same as accuracy
		indexLow = WHERE( spprec_out LT 1.0e-10, numLow )
		IF numLow GT 0 THEN spprec_out[indexLow] = !VALUES.F_NAN
		IF SIZE(spflux_out,/n_dim) EQ 2 THEN BEGIN
			spprecout = ( TOTAL( (spflux_out*spprec_out*good_err_mask) > 1e-12,2,/double, /NAN) / $
				TOTAL((spflux_out*good_err_mask) > 1e-12, 2, /double, /NAN) )

		ENDIF ELSE spprecout = (spflux_out*spprec_out*good_err_mask > 1e-12)
		; accuracy
		good_err_mask = (spflux_out GT 1e-12) ; byte mask

		; per EVE SDP meeting 4/20/11 Tom recommends a weighted average instead for estimating the accuracy
		indexLow = WHERE( spacc_out LT 1.0e-10, numLow )
		IF numLow GT 0 THEN spacc_out[indexLow] = !VALUES.F_NAN
		IF SIZE(spflux_out,/n_dim) EQ 2 THEN BEGIN
			spaccout = ( TOTAL( (spflux_out*spacc_out*good_err_mask) > 1e-12, 2, /double, /NAN) / $
				TOTAL((spflux_out*good_err_mask) > 1e-12, 2, /double, /NAN) )

			spflux = TOTAL(spflux_out,2,/double)
		ENDIF ELSE BEGIN
			; only one integration
			spaccout = (spflux_out*spacc_out*good_err_mask) > 1e-12
			spflux = spflux_out
		ENDELSE
		; not strictly enforce sanity
		spaccout = spaccout > spprecout
		spprecout = spprecout < spaccout
	ENDELSE

	;stop
	RETURN
END
