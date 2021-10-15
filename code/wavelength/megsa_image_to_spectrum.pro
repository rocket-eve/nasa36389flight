;+
; NAME:
;  megsa_image_to_spectrum
;
; PURPOSE:
;  Create a constant wavelength scale spectrum from a MEGS-A image
;
; CATEGORY:
;  Level 1 MEGS-A
;
; CALLING SEQUENCE:
;  megsa_image_to_spectrum, img, spwave, spflux,
;      mask=mask, err=err, good=good, imgwave=imgwave, sperr=sperr,
;     linefile='thelinefile'
;
; INPUTS:
;  img: a 2048x1024 MEGS-A image (ff & dark corrected)
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
;  megsa_image_to_spectrum_cal:
;   default_mask - the default mask from a saveset
;   default_imgwave - a default wavelength scale that should work for 
;                     center pointing
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;  megsa_image_to_spectrum, image, wave, flux
;   ; This example uses the default wavelength scale.
;
; MODIFICATION HISTORY:
;  03/22/07 DLW Original file creation
;
;-
pro megsa_image_to_spectrum, imgfull, spwave, spflux, $
  imgmask=imgmask, imgerr=imgerr, imgwgd=imgwgd, imgwave=imgwave, $
  sperr=sperr, linefile=linefile

; imgfull is a 2048x1024 ff corrected image from megs-b
; imgmask is a 2048x1024 mask to apply to the image (if undefined, a default mask is used)
; imgerr is a 2048x1024 float array of relative uncertainties to propagate
; imgwgd is a 2048x1024 bytarr of 0 (bad) or 1 (good) indicating where imgfull is good
; imgwave is an optional 2048x1024 wavelength image to use
;  (actual wavelengths are returned if not specified)
; spwave is the one-dimensional wavelength array used in regridding
; spflux is the output spectral flux
; sperr is the propagated uncertainty
; linefile=linefile allows the user to specify a file of lines to use for fitting (overriding any imgwave specified)

common megsa_image_to_spectrum_cal, default_mask, default_imgwave

if size(default_mask,/type) eq 0 then begin
    ;setup common block
;    restore,'megsa_default_mask.sav' ;restores default_mask variable
;    restore,'megsa_full_wave.sav' ;restores megsa_wave
    default_mask = intarr(2048,512)>1
    x=findgen(2048)
    wave1 = 3.83249 + (0.0135018 + x*(1.57669e-6))*x
    megsa_wave = rebin(wave1,2048,512)
    default_imgwave = temporary(megsa_wave) ;assign to default_imgwave variable
endif

; is a mask provided?
if size(imgmask,/type) eq 0 then imgmask = default_mask

; get the stripe image
;stripeimg = make_megsb_stripe( imgfull, imgmask )
stripeimg = imgfull[*,0:511]*imgmask

if size(imgwave,/type) eq 0 then begin
    ;imgwave is undefined, may need to do the full fit
    if size(linefile,/type) eq 0 then begin
        ; use default_imgwave
        print,'MEGSA_IMAGE_TO_SPECTRUM using default wavelengths'
        ;thisimgwave = make_megsa_stripe( default_imgwave, imgmask )
        thisimgwave = default_imgwave
        ; ERRORS ARE CURRENTLY IGNORED
        make_megsa_wave_img, stripeimg, thisimgwave, spwave, spflux ;, err=err
    endif else begin
        ;do the full fit with a specific linefile (the works!)
        print,'MEGSA_IMAGE_TO_SPECTRUM fitting wavelengths'
        make_megsa_wave_img, stripeimg, imgwave, spwave, spflux, file=linefile
    endelse
endif else begin
    ; get the stripe waveimg if needed
        print,'MEGSA_IMAGE_TO_SPECTRUM using wavelength argument (not a fit)'
    ;if n_elements(imgwave) ne n_elements(stripeimg) then $
      ;thisimgwave = make_megsb_stripe( imgwave, imgmask ) $
    ;else thisimgwave = imgwave

    if n_elements(imgwave) ne n_elements(stripeimg) then $
      thisimgwave = default_imgwave * imgmask $
    else thisimgwave = imgwave

    ; use the input wavelength scale instead of fitting
    n_spec = n_elements(stripeimg[0,*])
    n_wave = n_elements(stripeimg[*,0])

    if n_elements(imgwgd) ne n_elements(stripeimg) and size(imgwgd,/type) ne 0 then $
      thisimgwgd = imgwgd * imgmask $ ;make_megsa_stripe( imgwgd, imgmask ) $
    else thisimgwgd = bytarr(n_spec,n_wave)>1b

    if n_elements(imgerr) ne n_elements(stripeimg) and size(imgerr,/type) ne 0 then $
      thisimgerr = imgerr * imgmask $ ;make_megsb_stripe( imgerr, imgmask ) $
    else thisimgerr = fltarr(n_spec,n_wave)>1.

    ;spwave = dindgen(3100)*.025d0 + 31.d0
    spwave = dindgen(3600)*.01d0 + 3.d0
    n_out  = n_elements(spwave)
    spflux_out = dblarr(n_out)
    if size(sperr,/type) eq 0 then sperr=thisimgerr<1.
;stop
    regrid_spectra, thisimgwave, stripeimg, n_spec, n_wave, thisimgwgd, $
                    thisimgerr, fltarr(n_spec,n_wave)>1., $
                    spwave, spflux_out, n_out, wg_out, em_out, ea_out, $
                    sperr, fltarr(n_spec,n_wave)>1., $
                    acoro, dcoro, status

    spflux=total(spflux_out,2,/double) ; >.1

endelse

;stop
return
end
