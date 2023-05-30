pro make_megsb_wave_img, img, wimg, wout,spout, file=file

debug=0
;debug=1

if size(file,/type) eq 0 then file='megsb_solar_lines.dat'
;d=read_dat(file)
;keepers=where(d.cal ne 0,n_keepers)
;lines=d[keepers]

if size(wimg,/type) eq 0 then begin

    vctr=n_elements(img[0,*])/2
    wimg=(img<0)>0              ;for size

    striphgt=7
;throw out top and bottom to ensure lines are good on all strips
; use 100 pixels as the slit image height for the lines
; (less than actual, but need only the best part for fitting)
    slitimghgt=100
    n_strips=long(slitimghgt/striphgt)
    strip_sp=dblarr(2048,n_strips)
    strip_w = strip_sp

;loop over each strip to determine sp and wavelengths
; horizontal spectrum fit
    if debug ne 0 then begin
        !p.thick=1
        plot,[1,2],[1,2],xr=[36,38],xs=1,yr=[1e-5,1],/ylog,/nodata
        ;plot,[1,2],[1,2],xr=[73,75],xs=1,yr=[1e-5,1],/ylog,/nodata
    endif

    for i=0L,n_strips-1 do begin
        lo=(vctr-(slitimghgt/2)) + i*striphgt
        hi=lo+striphgt-1
        sp=total(img[*,lo:hi],2)
        w=megs_sp_wave_cal2(sp,err,s,file=file,/megsb,debug=0) ;use normal fitting for sp
        strip_sp[*,i]=sp
        strip_w[*,i]=w
        wimg[*,(lo+hi)/2]=w
        if debug ne 0 then $
          oplot,w,sp/max(sp),ps=10,co=i*10
    endfor
;stop
;loop over vertical columns to fit the wavelength scale between the strips
    vert=lindgen(n_elements(wimg[0,*]))
    for i=0,2047 do begin
        xpos=where(wimg[i,*] gt 0,n_xpos)
        coef=poly_fit(xpos,wimg[i,xpos],2)
        wimg[i,*]=coef[0] + coef[1]*vert + coef[2]*vert*vert; + coef[3]*vert*vert*vert
    endfor
endif

; determine a high-resolution spectrum (w and sp)
n_spec=n_elements(img[0,*])
n_wave=n_elements(img[*,0])
;wout=dindgen(3100)*.025d0 + 31.d0 ;min(wimg) 0.25 angstrom
wout=dindgen(3700)*.020d0 + 33.d0 ;min(wimg) 0.2 angstrom
n_out=n_elements(wout)
spout=dblarr(n_out)
;regrid_spectra, wimg, img>(1.e-14), n_spec, n_wave, bytarr(n_spec,n_wave)>1b, $
regrid_spectra, wimg, img, n_spec, n_wave, bytarr(n_spec,n_wave)>1b, $
                fltarr(n_spec,n_wave)>1., fltarr(n_spec,n_wave)>1., $
                wout, flux_out, n_out, wg_out, em_out, ea_out, $
                fltarr(n_spec,n_wave)>1., fltarr(n_spec,n_wave)>1., $
                acoro, dcoro, status
spout=total(flux_out,2,/double)>(1.e-9) ;.1
;stop
if debug eq 0 then return

;
; DEBUG ONLY
;


;!p.multi=[0,1,4]
lines=calc_megs_fwhm(wout,spout)
stop
sprel=spout/max(spout)
range=20 ;nm
for i=0,3 do begin
    lo=i*range + 30
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
