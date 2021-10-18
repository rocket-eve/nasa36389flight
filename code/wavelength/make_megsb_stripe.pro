;+
; NAME:
;  make_megsb_stripe.pro
;
; PURPOSE:
;  Shuffle the main diagonal of a MEGS-B image to a smaller rectangle image 
;  for use in determining a wavelength scale. Can also take a reduced stripe
;  image and expand it back to a full sized image. Off-diagonal data is zeroed.
;
; CATEGORY:
;  Level 1 MEGS-B
;
; CALLING SEQUENCE:
;  result_img = make_megsb_stripe( img [, mask] [,/undo] )
;
; INPUTS:
;  img: a 2048x1024 MEGS-B image (or /undo set, a 2048xN image to expand)
;
; OPTIONAL INPUTS:
;  mask: a 2048x1024 multiplicative mask of 1 along main image diagonal
;      : if mask is not present, a default mask will be used
;
; KEYWORD PARAMETERS:
;  /undo: set the undo switch to create a full size MEGs-B image from a 
;         stripe image (the reverse of the normal process)
;
; OUTPUTS:
;  result_img: the resulting image from applying the mask and shuffling
;            : if /undo is set, then result_img is a full-size image
;
; COMMON BLOCKS:
;  make_megsb_stripe_cal
;    default_mask: the mask restored from a saveset
;
; SIDE EFFECTS:
;  Process is not completely reversible, since off-diagonal data cannot be 
;  recovered.
;
; RESTRICTIONS:
;  When in doubt, use the default mask.
;
; PROCEDURE:
; 10) Setup common block data if necessary
; 20) Check parameters
; 30) Determine if this is an undo or not
;  40) apply the mask and define the output array to have the right type
;  50) loop over each column and copy the stripe of good data to out
;  60) resize output variable
; 70) return output variable
;
; EXAMPLE:
;  strip = make_megsb_stripe( img )
;  Then pass strip to make_megsb_wave_img.pro
;
; MODIFICATION HISTORY:
;  10/06/06 DLW Added comments.
;-

function make_megsb_stripe, in_img, mask, undo=undo

; reduce megsb image to a narrow stripe that contains only the slit images
; in_img is a 2d array (say 2048x1024)
; mask is the same dimension
; /undo option reverses the process to enlarge the shuffled image back 
;  to the size of the mask

common make_megsb_stripe_cal, default_mask

;
; 10) Setup common block data if necessary
;
if size(default_mask,/type) eq 0 then begin
    local_path=file_dirname(file_which('megsb_image_to_spectrum.pro'))+'/'
    restore,local_path+'megsb_default_mask.sav' ;default_mask is a named variable
    ; default mask was pulled from Neon spectrum analysis
endif

;
; 20) Check parameters
;
if n_params() gt 2 or n_params() lt 1 then begin
    print,'USAGE: stripe = make_megsb_stripe( img [,mask] [/undo] )'
    print,' img is a 2048x1024 image'
    print,'  stripe is a 2048x200 image (unless a differen mask is applied)'
    print,'  /undo will reverse the process'
    return,-1
endif

if n_params() eq 1 then mask = default_mask ;else mask was specified

;
; 30) Determine if this is an undo or not
;
if keyword_set(undo) eq 0 then begin
    ;
    ; keyword was not set
    ;
    ; 40) apply the mask and define the output array to have the right type
    ;
    img = in_img * mask
    out = make_array(size(img,/dim),type=size(img,/type)) ; associate with the same data type
    maxdist=0
    ;
    ; 50) loop over each column and copy the stripe of good data to out
    ;
    for i=0L,n_elements(img[*,0])-1 do begin
        good=where(img[i,*] gt 0, n_good)
        if n_good gt 1 then begin
            low=min(good,max=hi)
            dist=hi-low
            out[i,0:dist]=img[i,low:hi] ;should be a simple translation
            maxdist = maxdist>dist
        endif
    endfor
    ;
    ; 60) resize output variable
    ;
    out=out[*,0:maxdist]
endif else begin
    ; undo the above shuffling
    out = make_array(size(mask,/dim),type=size(in_img,/type)) ; define same type
    for i=0L,n_elements(out[*,0])-1 do begin
        good=where(mask[i,*] gt 0,n_good)
        if n_good gt 1 then begin
            low = min(good,max=hi)
            dist = (hi-low) < (n_elements(in_img[0,*])-1)
            ;out[i,low:hi]=in_img[i,0:dist]
            out[i,low:(low+dist)]=in_img[i,0:dist]
        endif
    endfor
endelse

;
; 70) return output variable
;
return,out
end
