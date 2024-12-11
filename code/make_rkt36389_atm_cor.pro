;+
;-
pro make_rkt36389_atm_cor

numberstr='389'
datestr='May 3, 2023'
yd = 2023123
; launch at 18:30 UT (reltime=0)
sodoffset = 18.d0*60.*60. + 60*30. ; seconds of day at reltime=0

;raw=read_dat('36'+numberstr+'_radar_data.txt')
;raw = read_dat('RADAR_36.336_06182018.dat') ; ~50 ms sampling
;raw = read_dat('Woods_36-353_Radar.dat') ; ~50 ms sampling (TEMPORARY)
;raw = read_dat('Radar_WOODS_36.353_PAYLOAD_09092021.dat') ; ~50 ms sampling
raw = read_dat('Radar_Data_36_389.dat') ; ~50 ms sampling
orig=raw
;
; Header
;

;                 TIME       SLT RG(M)    AZIMUTH    ELEVATION   HOR RG(M)   NORTH SOUTH EAST WEST   V. AZ     ALT.(M)      LATITUDE    LONGITUDE      TVEL(M).    E/W VEL.    N/S VEL.    ALT VEL(M)   FLT EL     FLT AZ      LK RG(M)    LOOK AZ     LOOK EL
;
;15022  lines:  36.353 Radar processed data for Payload
;21 columns: see header for definitions:  time=data[1,*], range_m=data[4,*], alt_m=data[9,*]

;:  NASA 36.353 Radar data for Sept 9 launch at 17:25:00 UT

; 36.353 file is 50 ms (units sec)
; 36.336 file is 100 millisecond sampling initially, then larger jumps
; 36.290 file is 50 millisecond sampling
; interpolate raw to a uniform sampling of 0.1 seconds

maxaltm = max(raw[9,*],wmax)
print,''
print,'NASA 36.'+numberstr+' apogee of '+strtrim(maxaltm/1000.,2)+' km at T+'+strtrim(raw[1,wmax],2)+' seconds'
print,''


; change time from ms to seconds
raw = 1.d * (raw) ; enforce double precision
time_sec = raw[1,*] ;double(reform(raw[0,*])) / 1000.d ; ms to sec
raw[0,*] = time_sec ; replace nonused column 0 with something usefull
;altitude_m = raw[9,*]
;raw[9,*] = raw[9,*] / 100.d ; cm to m

;n_samples = max(raw[1,*]) ; time span from 0 in seconds
n_samples = max(time_sec) ; time span from 0 in seconds
;n_samples = long(n_samples*10.d) ; 10 samples per second
n_samples = long(n_samples) ; 1 sample per second
n_cols = n_elements(raw[*,0])
newraw = dblarr( n_cols, n_samples )

;newraw[1,*] = raw[1,0] + dindgen(n_samples)/10.d ; 0.1 sec resampling
;newraw[1,*] = raw[1,0] + dindgen(n_samples)/1.d ; 1.0 sec resampling
newraw[0,*] = time_sec[0] + dindgen(n_samples)/1.d ; 1.0 sec resampling
newraw[1,*] = newraw[0,*] ; copy time

;for col=2L,n_cols-1 do newraw[col,*] = interpol(raw[col,*],raw[1,*], newraw[1,*])
for col=2L,n_cols-1 do newraw[col,*] = interpol(raw[col,*], time_sec, newraw[1,*])
oldraw = temporary(raw)
raw = temporary(newraw); uniform time sampling
time_sec = raw[1,*] ; choose same spacing as raw
stop


;;skipsize = round( 10./(raw[1,2]-raw[1,1]) ) ; 1 sample per 10 seconds
;skipsize = round( 10./(time_sec[2]-time_sec[1]) ) ; 1 sample per 10 seconds

;apogee
;apogeekm=max(raw[9,*],wmax)
;apogeereltime=raw[1,wmax]-raw[1,0]
apogeekm=max(raw[9,*],wmax)
apogeereltime = time_sec[wmax] - time_sec[0]
print,'***'
print,'NASA 36.'+numberstr+' apogee of '+strtrim(apogeekm/1000.,2)+' km at T+'+strtrim(apogeereltime,2)+' seconds'
print,'***'
; for 36.258 the first integration is at 1.0 seconds, so chop to 10-sec samples

junk=min(abs(time_sec - 8.2),idx) ; start at 8.2 sec and skip by 10 sec
;filtered=raw[*,idx:*:200] ; chop off data before first CCD integration in image list
filtered=raw[*,idx:*:10] ; chop off data before first CCD integration in image list
;ft7a=76.6
;ft7=81.6
;fap=24.0
ft7 = 158.7 ; adj
ft7a = 111.4 ;152.0 ; 3/24 through 5/23
; download date range of data with lisird and run  more penticton_radio_flux.csv | tail -183 | awk -F',' '{print $2}' | awk -f ~/bin/avg.awk
fap = 3.0 ; Potsdam https://kp.gfz-potsdam.de/en/data

;YYYY MM DD Days  Days    Bart Db  Kp0    Kp1    Kp2    Kp3    Kp4    Kp5    Kp6    Kp7    ap0  ap1  ap2  ap3  ap4  ap5  ap6  ap7    Ap  SN   F10obs   F10adj D
;2023 05 03 33360 33360.5 2587 25  0.000  0.000  0.333  0.333  1.000  1.000  0.667  2.000    0    0    2    2    4    4    3    7     3 124    156.2    158.7 2


altm=reform(filtered[9,*]) ; altitude in m above sea level
altkm=altm/1000. ; m to km
reltime=reform(filtered[0,*]) ; seconds since launch
latdeg=reform(filtered[10,*]) ; around 32-33 deg
londeg=reform(filtered[11,*]) ; around -106.3 deg
;latdeg=reform(filtered[7,*]) ; around 32-33 deg
;londeg=reform(filtered[8,*]) ; around -106.3 deg
;lookaz=reform(filtered[19,*])
;lookel=reform(filtered[20,*])

apogeekm = max(altkm,apogeeidx)

; call SEE code to get atmospheric absorption for wavelengths we care
; about
wave = findgen(2000)*0.1 + 0.05 ; (0.1 nm bins)

; what is solar zenith angle?
;  https://www.esrl.noaa.gov/gmd/grad/solcalc/azel.html
; Launch pad is at 32.41782, -106.319702
; during launch, stays between 32.41 and 33.1, and -106.3 to -109.146
; (almost linear change)
; indicates solar elevation is 72.04 at 18:32UTC so sza is 90-72.04=17.96
 ; at 18:32 UTC 32deg 25' 4'', +106deg 18' 10'' (radar file)
; solar elev is 73.04deg at 18:47 UT, sza=90-73.04=16.96 deg
; I assume that 17.96 deg is pretty close for all solar observations
;sza = 17.5 ; mean of 17 and 18 ; 36.258

; for 36.389 May 3, 2023
; solar elevation is 71.76 at 18:35, so sza=90-71.76=18.24 deg
sza=18.24


; for 2021 sep 9, solar elevation is 58.95 deg
; so sza is 90-58.95=31.05
;sza=31.05

;;sza= 9.12
;; for Oct 2013, sun is low in the sky
;; used 32deg 25' North and 106deg 14' West (positive on page)
;; set UTC offset=0 and use UTC 18:00:00 and 18:15:00
;; solar elevation is 45.07 and 45.87, so sza is 90-45.07 = 44.93 or 44.13
;sza = 44.53 ; mean of 44.93 and 44.13 near center time

datarec={reltime:0.d, altkm:0., latdeg:0., londeg:0., trans:dblarr(n_elements(wave))}
result={rkt:'36.'+numberstr, ft7a:ft7a, ft7:ft7, fap:fap, sza:sza, wave:wave, data:replicate(datarec,n_elements(reltime))}

result.data.reltime=reltime
result.data.altkm=altkm
result.data.latdeg=latdeg
result.data.londeg=londeg

trans = fltarr(n_elements(reltime),n_elements(wave))

;!p.multi=[0,1,2]
!p.multi=0
!p.charsize=2
;!p.background=0
;!p.color='ffffff'x
!p.background='ffffff'x & !p.color=0
loadct,39
device,decomp=0
plot,reltime,altkm,xtit='Seconds since Launch',ytit='Altitude (km)',tit='NASA 36.'+numberstr+' '+datestr,xs=1,xr=[0,600],ps=-4
oplot,apogeekm+[0,0],!y.crange,lines=1
;oplot,212+[0,0],!y.crange,lines=1
;oplot,332+[0,0],!y.crange,lines=1

altthick=20
for i=150L,300,altthick do begin
  gd = where( altkm gt i and altkm le i+altthick, n_gd )
  if n_gd gt 0 then begin
     thisco = (i-150)*255./(300-150.)
     oplot,reltime[gd],altkm[gd],ps=1,symsize=2,co=thisco
     for j=0,n_gd-1 do oplot,reltime[gd[j]]*[1,1], [0,altkm[gd[j]]],co=thisco,thick=2
  endif
endfor
stop

firstplot=0

print,'Apogee conditions'
print,' sec of day       : ',sodoffset + apogeereltime
print,' relative seconds : ',apogeereltime
print,' latitude  (deg): ',latdeg[apogeeidx]
print,' longitude (deg): ',londeg[apogeeidx]
print,' altitude (km)  : ',altkm[apogeeidx]
print,' F10.7   : ',ft7
print,' <F10.7> : ',ft7a
print,' Ap      : ',fap
print,' SZA     : ',sza

;
; pin MSIS result to match ESP 30 channel at sodoffset+apogeereltime
; Andrew says ESP absorption is 4.1% at 30 at apogee
; explore msis absorption at 30.378 by changing f10.7 and <f10.7>
; there will be a curve of solutions, prefer changing <f10.7>
;
; w30 is an array because get_atmos_corr needs wavelength to be an array
w30 = 30.378 + (.01*[-1.,0.,1.]) ; nm 3 element array
; all atomic cross sections are smooth, to w30 does not ahve to be exact
; and atomic oxygen is the most important constituent near apogee
espPercent = 1.0 - 0.041 ; fraction from 1.
halfwidth = 60
f10range  = ft7  + (findgen(halfwidth*2) - halfwidth)
f10arange = ft7a + (findgen(halfwidth*2) - halfwidth)
t30 = fltarr(n_elements(f10range), n_elements(f10arange))
for i=0,n_elements(f10range)-1 do begin
   for j=0,n_elements(f10arange)-1 do begin
      ; get_atmos_corr needs wavelength to be an array
      get_atmos_corr, yd, sodoffset + apogeereltime, sza, latdeg[apogeeidx], londeg[apogeeidx], altkm[apogeeidx], f10arange[j], f10range[i], fap, w30, t30tmp, err, status
      t30[i,j] = t30tmp[1]
   endfor
endfor

!p.multi=0
!p.background='ffffff'x & !p.color=0
levels = findgen(10)*0.01-0.05 + espPercent
annot = strtrim(strcompress(string(levels,form='(f6.3)')))
contour, t30, f10range, f10arange,xs=1,ys=1, c_charsize=2, $
         xtit='F10',ytit='<F10.7>',tit='NRL-MSISe00 at 30.378 nm, Apogee='+strtrim(string(apogeekm,form='(f7.1)'),2)+', AP='+strtrim(string(fap,form='(f5.1)'),2)+' F10='+strtrim(string(ft7,form='(f5.1)'),2)+' <F10>='+strtrim(string(ft7a,form='(f5.1)'),2)+' Target=0.959', $
         c_annotation = annot, levels = levels
oplot,[1,1]*ft7,[1,1]*ft7a,ps=1,co='fe'x,symsize=3

; interpolate along F10 to find <F10> that produces the ESP measurement
newft7a = interpol(f10arange, reform(t30[halfwidth,*]), espPercent)
oplot,[1,1]*ft7, [1,1]*newft7a,co='cc00'x,ps=2,symsize=3
print,' interpolated <F10.7>=',newft7a
stop

print,'replacing old ft7a with newft7a to match ESP at 30nm'

ft7a = newft7a

stop
;


for i=0L,n_elements(altkm)-1 do begin
   ; 17:25 
  ;get_atmos_corr, yd, 17.d0*60.*60. + 60*25. + reltime[i], sza, latdeg[i], londeg[i], altkm[i], ft7a, ft7, fap, wave, transmission, err, status

   ; 18:30
   get_atmos_corr, yd, sodoffset + reltime[i], sza, latdeg[i], londeg[i], altkm[i], ft7a, ft7, fap, wave, transmission, err, status

  trans[i,*]=transmission
  result.data[i].trans = transmission

  if firstplot eq 0 then begin
     plot,wave,transmission,yr=[0.6,1],ys=1,xs=1,xr=[0,110],xtit='Wavelength (nm)',ytit='Transmission',tit='MSIS00 Transmission at all altitudes',/nodata
     firstplot=1
  endif
  if altkm[i] ge 100 then begin
    oplot,wave,transmission,color=(((altkm[i]-150.)>0) * 255.)/(300-150.)
    ;stop
  endif
endfor

outfile='rocket_36'+numberstr+'_atmtrans.ncdf'
write_netcdf,result,outfile,s,/clobber
print,'wrote '+outfile
stop

!p.multi=[0,1,2]
!p.charsize=2
!p.background=0
!p.color='ffffff'x
loadct,39
device,decomp=0
plot,reltime,altkm,xtit='Seconds since Launch',ytit='Altitude (km)',tit='NASA 36.'+numberstr+' '+datestr,xs=1,xr=[100,450],ps=-4
oplot,apogeekm+[0,0],!y.crange,lines=1

altthick=20
for i=150L,300,altthick do begin
  gd = where( altkm gt i and altkm le i+altthick, n_gd )
  if n_gd gt 0 then begin
     thisco = (i-150)*255./(300-150.)
     oplot,reltime[gd],altkm[gd],ps=1,symsize=2,co=thisco
     for j=0,n_gd-1 do oplot,reltime[gd[j]]*[1,1], [0,altkm[gd[j]]],co=thisco,thick=2
  endif
endfor
for i=0L,n_elements(altkm)-1 do begin

  get_atmos_corr, yd, sodoffset + reltime[i], sza, latdeg[i], londeg[i], altkm[i], ft7a, ft7, fap, wave, transmission, err, status

  trans[i,*]=transmission
  result.data[i].trans = transmission

  if i ge 13 and i le 35 then begin
    if i eq 13 then begin
      plot,wave,transmission,yr=[0.,1],ys=1,xs=1,xr=[0,100],xtit='Wavelength (nm)',ytit='Transmission',tit='MSIS00 Transmission at all altitudes'
    endif
    oplot,wave,transmission,color=(((altkm[i]-150.)>0) * 255.)/(300-150.)
    ;stop
  endif
endfor
stop

return
end
