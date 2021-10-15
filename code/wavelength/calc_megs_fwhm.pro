function calc_megs_fwhm, w, sp ;, width=width

; pass in wavelength scale in nm, and spectrum
; returns a structure with a record for each peak
;  contains - w: wavelength center of peak
;           - pixel: pixel fraction of center
;           - max:  peak value
;           - width: peak width at half max (in nm)

linerec={w:0.0d0, pixel:0.0d, max:0.0d0, width:0.0d0}

if max(sp) lt 2 then spectrum=sp*1.d7 else spectrum=sp
if max(sp) lt 1000 then spectrum=sp*1.d3 else spectrum=sp

;call Tom's routine to extract lines
list=extract_megs_lines(spectrum) ;,width=width)
;stop
lines=replicate(linerec,n_elements(reform(list[0,*])))
lines.pixel = reform(list[0,*])
if max(sp) lt 2 then $
  lines.max = reform(list[1,*])/1.d7 else $
    lines.max = reform(list[1,*])
lines.w     = interpol(w,dindgen(n_elements(w)),reform(list[0,*]))
; width is in pixels
lines.width = reform(list[2,*])
; put onto nm scale
for i=0L, n_elements(list[0,*]) - 1 do begin
    tmp = min( abs(lines[i].w - w), idx ) ;get closest index in w for peak
    delta = w[idx+1] - w[idx] ; pixel spacing in nm
    lines[i].width = lines[i].width * delta
endfor
;stop
return,lines
end
