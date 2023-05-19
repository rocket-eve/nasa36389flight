pro make_megsa_wave_img, img, wimg, wout,spout, file=file, slit1=slit1

; break image into horizontal strips of spectra
; fit each spectrum (strip) individually, 
; vertically fit the fit (strip) results
;  and interpolate/extrapolate to the rest of the columns
; use regrid spectra to get a fixed wavelength
; scale for all spectra

debug=1
;if keyword_set(slit1) then debug=1

if size(file,/type) eq 0 then file='megsa_solar_lines.dat'

if size(wimg,/type) eq 0 then begin

    ;vctr=n_elements(img[0,*])/2
    ; optimize vctr to center of mass of horizontal sum
    t=total(img,1)
    t=t-min(t)
    vctr=round(total(t*findgen(n_elements(img[0,*]))) / total(t))
    ; or just use 299 for center pointing
    print,'vctr='+strtrim(vctr,2)
    wimg=img*0.              ;for size

    if keyword_set(slit1) then $
       striphgt=7 $ ;15  $           ;7 ;30 ;15
    else striphgt=7
;throw out top and bottom to ensure lines are good on all strips
; use 150 pixels around CM as the slit image height for the lines
; (less than most lines, but we want to use only the best part for fitting)
    if keyword_set(slit1) then $
       slitimghgt=150 $           ;100 ;150
    else slitimghgt = 100
    n_strips=long(slitimghgt/striphgt)
    strip_sp=dblarr(2048,n_strips)
    strip_w = strip_sp

;loop over each strip to determine sp and wavelengths
; horizontal spectrum fit
    if debug ne 0 then begin
        !p.thick=1
    ;plot,[1,2],[1,2],xr=[102,104],xs=1,yr=[1e-5,1],/ylog,/nodata
        plot,[1,2],[1,2],xr=[0,40],yr=[1e-5,1],/ylog,/nodata
    endif

    for i=0L,n_strips-1 do begin
        lo=(vctr-(slitimghgt/2)) + i*striphgt
        hi=lo+striphgt-1
        sp=total(img[*,lo:hi],2)
;        w=megs_sp_wave_cal(sp,err,s,file=file,debug=0,slit1=slit1) ;use normal fitting for sp
        w=megs_sp_wave_cal2(sp,err,s,file=file,debug=0,slit1=slit1) ;use normal fitting for sp
        strip_sp[*,i]=sp
        strip_w[*,i]=w
        wimg[*,(lo+hi)/2]=w
        if debug ne 0 then begin
           oplot,w,sp/max(sp),ps=10,co=i*255./n_strips
           print,lo,hi
        endif
        ;if debug ne 0 then stop
    endfor
;stop
;loop over vertical columns to fit the wavelength scale between the strips
    vert=lindgen(n_elements(wimg[0,*]))
    for i=0,2047 do begin
        xpos=where(wimg[i,*] gt 0,n_xpos)
        ; quadratic
        ;coef=poly_fit(xpos,wimg[i,xpos],2)
        ;wimg[i,*]=coef[0] + coef[1]*vert + coef[2]*vert*vert

        ; linear 
        coef=poly_fit(xpos,wimg[i,xpos],1)
        wimg[i,*]=coef[0] + coef[1]*vert ; + coef[2]*vert*vert
    endfor
endif

; determine a high-resolution spectrum (w and sp)
n_spec=n_elements(img[0,*])
n_wave=n_elements(img[*,0])

wout=dindgen(3600)*.010d0 + 3.d0 ;min(wimg)
n_out=n_elements(wout)
spout=dblarr(n_out)
regrid_spectra, wimg, img, n_spec, n_wave, bytarr(n_spec,n_wave)>1b, $
                fltarr(n_spec,n_wave)>1., fltarr(n_spec,n_wave)>1., $
                wout, flux_out, n_out, wg_out, em_out, ea_out, $
                fltarr(n_spec,n_wave)>1., fltarr(n_spec,n_wave)>1., $
                acoro, dcoro, status
spout=total(flux_out,2,/double)>(1.e-9) ;.1

;stop
return

if debug eq 0 then return

;
; DEBUG ONLY
;


;!p.multi=[0,1,4]
stop
lines=calc_megs_fwhm(wout,spout)
stop
sprel=spout/max(spout)
range=10 ;nm
for i=0,3 do begin
    lo=i*range + 10
    hi=lo+range
    tmp=where(wout ge lo and wout le hi,n_tmp)
    if n_tmp gt 0 then begin
        plot,wout,sprel,xr=[lo,hi],xs=1, ps=10, $
             yr=[min(sprel[tmp])>(1e-4),max(sprel[tmp])*2], $
             /ylog,xtit='Wavelength (nm)',ytit='relative signal'
        for j=0,n_elements(lines)-1 do begin
            oplot,[1,1]*lines[j].w,10^!y.crange,lines=1
            xyouts,lines[j].w,lines[j].max*1.04/max(spout),string(lines[j].w,form='(f7.3)'),orient=90.,charsize=2
        endfor
    endif
    stop
endfor

stop
return
end
