function remove_atm_absorption_36389, sptime, wave, sp_in

; sptime is a 1d array of times relative to launch
; sp is a 2-d array of [time, irradiance]
;  ONLY pass in the integration that are solar measurements!!!

@config36389
  
sp = sp_in ; working copy

dims=size(sp,/dim) ; get dimension lengths
if dims[0] ne n_elements(sptime) then begin ; does 0th dim length match?
  sp = transpose(sp) ; swap dimensions
  dims = size(sp,/dim)
  if dims[0] ne n_elements(sptime) then begin
    print,'ERROR: remove_atm_absorption - sptime and sp do not agree in length'
    help, sp, sptime
    stop
  endif
endif
ntime=dims[0]
nwave=dims[1]

routinename = ((scope_traceback(/str))[-1]).routine

; use data from msis00e run on timed-see
;print,'INFO: '+routinename+' using final 36.353 flight profile'
print,'INFO: '+routinename+' using final 36.389 flight profile'
workingdir = file_dirname(routine_filepath()) ; in code/compare_merge
datadir = file_dirname(file_dirname(workingdir))+'/data'
;read_netcdf,datadir+'/rocket_36353_atmtrans.ncdf',atm
read_netcdf,datadir+'/rocket_36389_atmtrans.ncdf',atm

!p.background='ffffff'x
!P.color=0

co=['black','midnight blue','indigo','purple','blue violet','blue','royal blue','teal','dark green','olive drab','medium sea green','yellow green','lime green','goldenrod','dark orange','coral','salmon','indian red','firebrick','crimson','deep pink','red','dark red','brown','chocolate','peru','tan','dark khaki','khaki','navajo white','blanched almond','light gray','silver','dark gray','gray','dim gray','light slate gray','dark slate gray']

;plot,atm.data.reltime,atm.data.altkm,xtitle='Time since launch (seconds)',$
;  ytitle='Altitude (km)',ps=-4,xr=[0,500],xs=1
;write_pngfile,'rktaltitudevstime.png' ; does not work anymore

; old black and white plot using just radar
;p=plot(atm.data.reltime,atm.data.altkm,xtitle='Time since Launch (seconds)', $
;  ytitle='Altitude (km)', title='NASA Rocket 36.336',$
;  xrange=[0,500], dimensions=[800,600],$
;  lines=1, symbol='circle' )
;for i=0L,ntime-1 do p=plot([1,1]*sptime[i],[0,1]*atm.data.altkm[i],color=co[i],/overplot,linestyle='solid')
;p.save, 'rktaltituderadarvstime.png',width=800,height=600
;stop
;p.close


newatmrec={reltime:0.d, altkm:0.d, trans:dblarr(n_elements(atm.data.trans[*,0]))}
newatm = replicate(newatmrec,ntime)

timestampoffset=5.0 ; if time is start, then add 5 sec to get center
;timestampoffset=-5.0 ; if time is end, then subtract 5 sec to get center

for i=0L,ntime-1 do begin
  for w=0L,n_elements(atm.wave)-1 do begin
    newatm[i].trans[w]=interpol( $
      atm.data.trans[w,*], atm.data.reltime, sptime[i]+timestampoffset, /lsquad)
 endfor
  newatm[i].reltime = sptime[i]
  newatm[i].altkm = interpol(atm.data.altkm, atm.data.reltime, sptime[i]+timestampoffset, /lsquad)
endfor

; need to divide the spectra by the transmission to "remove" the
; atmosphere
newsp = sp
allcor = fltarr(ntime, nwave)

for itime=0L,ntime-1 do begin
  ; find closest time in atm.data.reltime to sptime[itime]
  junk = min(abs(sptime[itime]-newatm.reltime), idx)
  ; use idx
  cor = interpol(reform(newatm[idx].trans), reform(atm.wave), wave, /lsquad)
  allcor[itime,*] = cor ; keep all corrections to make a figure
  newsp[itime,*] /= reform(cor)
endfor

; altitude plot
p=plot(newatm.reltime,newatm.altkm,xtitle='Time since Launch (seconds)', $
  ytitle='Altitude (km)', title='NASA Rocket 36.'+numberstr,$
  xrange=[min(sptime)-20,max(sptime)+20], dimensions=[800,600],$
  linestyle='solid', /buffer )
for i=0L,ntime-1 do begin
  p=plot([1,1]*sptime[i],[0,1]*newatm[i].altkm,color=co[i],/overplot,linestyle='solid',thick=10)
  p=plot([1,1]*sptime[i]+2,[0,1]*newatm[i].altkm,color=co[i],/overplot,linestyle='solid',thick=10)
  p=plot([1,1]*sptime[i]-2,[0,1]*newatm[i].altkm,color=co[i],/overplot,linestyle='solid',thick=10)
endfor
outpngfile='rktaltitudevstime.png'
p.save, outpngfile,width=800,height=600
print,'saved '+outpngfile
p.close


if max(wave) lt 40 then begin
  det='MEGS-A'
  name='megsa'
  yrange=[0.8,1]
endif else begin
  det='MEGS-B'
  yrange=[0.5,1]
  name='megsb'
endelse

p = plot( wave, allcor[0,*], $
  xtitle='Wavelength (nm)', ytitle='Absorption', $
  title='NRL-MSIS00E Transmission for '+det+' Rocket 36.'+numberstr, $
  xrange=[min(wave),max(wave)], yrange=yrange, /stairstep, symbol='none', $
  dimensions=[800,600], /buffer )

for i=1L,ntime-1 do p = plot( wave, allcor[i,*], /overplot, /stairstep, symbol='none', color=co[i] )

outpngfile='rktatmcorvswave_'+name+'.png'
p.save, outpngfile,width=800,height=600
print,'saved '+outpngfile
p.close

tmp = max(newatm.altkm, apogeeidx) ; need index to apogee

if name eq 'megsa' then begin
  ; plot 30.4 line total vs time
  lo=2725 & hi=2750
  sig = total( sp[*,lo:hi],2) ; 30.4 line signal in each image
  ; is this MEGS-A slit 2?
  if mean(sig) gt .03 then begin
    newsig = total(newsp[*,lo:hi],2)
    scale=mean(newsig[apogeeidx-1:apogeeidx+1])
    rsig=sig/scale
    rnewsig=newsig/scale
    xrange=[min(sptime)-20,max(sptime)+20]

    p = plot(sptime, rsig, yr=[.8,1.0+3.*stddev(rnewsig)],$
      xtitle='Time since launch (seconds)',$
      ytitle='Normalized irradiance',$
      title='36.'+numberstr+' 30.4 nm Measurement and with MSIS Correction',$
      xrange=xrange, dimensions=[800,600],$
      linestyle='solid', symbol='circle',/buffer )
    p = plot( xrange, [1,1], linestyle='dash', /overplot )
    p = plot( sptime, rnewsig, color='dark green', linestyle='solid', symbol='square', /overplot )
    for i=0L,ntime-1 do begin
      p=plot([1,1]*sptime[i],[0,1]*rsig[i],color=co[i],/overplot,linestyle='solid',thick=10)
      p=plot([1,1]*sptime[i]+2,[0,1]*rsig[i],color=co[i],/overplot,linestyle='solid',thick=10)
      p=plot([1,1]*sptime[i]-2,[0,1]*rsig[i],color=co[i],/overplot,linestyle='solid',thick=10)
    endfor
    p = text(/norm,.15,.6,'Measurement')
    p = text(/norm,.15,.78,'Corrected, 1-sigma spread='+strtrim(string(100.*stddev(rnewsig),form='(f6.2)'),2)+'%',color='dark green')
;stop
    outpngfile='rktatmcor304vstime_'+name+'.png'
    p.save, outpngfile,width=800,height=600
    print,'saved '+outpngfile
    ;stop
    p.close
 endif ; slit 2
endif else begin
  ; megsb

endelse

;stop
return,newsp
end


function unused_remove_atm_absorption


; alternative method, needs pretty clean data to derive an exponential
; fit, not really feasible for most wavelengths




; altkm is a 1d array with length of time
; perigee is the index into altkm
; perideeidx is the sp index into the time dimension used as the perigee

; use mpfit to determine the exponential

;d=read_dat('data/RADAR_36.336_06182018.dat')
radarfile=file_search('data/*RADAR*36.'+numberstr+'*',count=count)
if count eq 1 then d=read_dat(radarfile) else begin
   print,'ERROR: remove_atm_absorption_36336 - cannot locate radar data file'
   stop
endelse

altkm=d[9,*] / 1000.            ; file is m, convert to km
altsec = d[1,*] ; - d[1,0] ; sec since launch
maxalt=max(altkm, perigee)
junk = min(abs( altsec[perigee] - sptime), perigeeidx )

relsp = sp ; make a normalized copy of the spectra

; use transpose if dimensions are not right

dims=size(relsp,/dim)
if dims[0] ne n_elements(sptime) then begin
  relsp = transpose(relsp)
  dims = size(relsp,/dim)
  if dims[0] ne n_elements(sptime) then begin
    print,'ERROR: remove_atm_absorption - sptime and sp do not agree in length'
    stop
  endif
endif
ntime=dims[0]
nwave=dims[1]

; normalize to 1 at perigee
norm = mean(relsp[[perigeeidx-2:perigeeidx+2],*],dim=1)
for itime=0L,ntime-1 do begin
  relsp[itime,*] /= norm
endfor

; now call the fitter for each wavelength on the way up to perigee
; apply the correction using altitude for the descent
correction = dblarr(dims) + 1.d0 ; pre-populate with 1.0
for iwave = 0L, nwave-1 do begin
   ts = reform(relsp[0:perigeeidx,iwave])
   gd=where(ts gt 1e-6,n_gd)
   if n_gd gt 5 then begin
      spread = stddev(ts[gd])
      ; using the stddev is probably not the best thing
      if spread gt 1e-5 and spread lt .3 then begin
        plot,gd,ts[gd], tit=strtrim(iwave,2),yr=[0,1.1],ps=-4
        ; fit to y=1-exp(p0*(p1-x))
        stop
      endif
   endif
endfor

; corrected spectrum has the atmosphere removed
result = sp * correction
stop
return, result
end
