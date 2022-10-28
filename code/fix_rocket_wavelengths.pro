function mygaussnterms4, x, p
   ; match gaussfit nterms=4
  z = (x - P[1]) / P[2]
  return, P[0] * exp( (-0.5d0) * z * z ) + P[3]
end

function get_ref_lines

  ; probably need a reference set of solar lines to pin to absolute scale
  reflines = [33.27900d0, $         ; Al X -.06bins
              33.41783, $         ; Fe XIV +0.06
              33.54100, $         ; Fe XVI -0.22 bins
              ;34.19511, $         ; Si IX (Fe XVIII blend with S VII 34.2026) +0.38 suspicious
              34.2026, $ ; NIST S VII +0.002
              34.74026, $       ; Si X -0.05
              34.97917, $         ; Si IX 0.003
              ;;34.98602, $         ; Si IX
              35.60381, $         ; Si X 0.01
              36.07590, $         ; Fe XVI 0.03
              36.44670, $         ; Fe XII -0.01
              ;36.45039, $         ; Si XI
              36.80713, $         ; Mg IX (blend with Fe V 36.7992) -0.28
              37.40040, $         ; O III -0.02
              ;37.40730, $         ; O III -0.37
              ;37.93080, $         ; Ne III
              38.40310, $         ; C IV (Fe XX) -0.01
              ;38.73560, $         ; N IV
              ;38.91360, $         ; Ar XVI
              ;39.24330, $         ; Al IX 0.65 suspicious
              ;39.55570, $       ; O III 
              ;39.5633, $        ; Fe V NIST +0.08
              ;40.19280, $         ; Ne VI -0.48 too blended
              40.32700, $         ; Ne VI -0.13
              41.11660, $         ; Na VIII (Fe XVIII) +0.13
              41.72580, $         ; Fe XV 
              ;41.76611, $       ; S XIV 7.1
              41.97150, $         ; C IV -0.19
              43.04550, $         ; Mg VIII +0.23
              43.31730, $         ; Ne VI -0.17
              43.49239, $         ; Mg VII (O III) +0.04
              43.67340, $         ; Mg VIII -0.15
              ;43.91771, $         ; Mg IX -0.81
              43.9048, $        ; O IV 
              ;44.12000, $         ; Mg IX maybe?
              ;44.39737, $         ; Mg IX -0.32 too blended
              ;44.57011, $         ; S XIV 7.1 
              45.96270, $         ; C III -0.17
              46.52200, $         ; Ne VII -0.03
              46.98250, $         ; Ne IV 0.028
              ;48.13740, $         ; Ne V -0.22
              48.13630, $         ; Ne V -0.16
              48.29970, $         ; Ne V -0.11
              ;48.94950, $         ; Ne III 0.24 blend with 48.9629 0.24
              48.9562, $ ; BLEND Ne III (48.9495+48.9629)/2.
              49.94066, $         ; Si XII -0.03
              ;50.81780, $         ; O III -1.14 blend with 50.768
              50.79290, $         ; BLEND O III (50.8178+50.7680)/2. 0.10
              50.997, $           ; NIST He I +0.12
              ;51.00429, $       ; Fe XIII -0.24
              ;51.708, $         ; NIST He I Not detectable
              ;51.771, $         ; NIST He I Too dim
              ;51.863, $         ; NIST He I Too dim
              51.207, $         ; NIST He I +0.006
              51.5596, $         ; NIST He I +0.04
              ;51.56180, $         ; CHIANTI He I
              52.06661, $         ; Si XII -0.03
              ;52.22140, $         ; He I
              52.2186, $         ; NIST He I -0.14
              52.57940, $         ; O III 0.02
              ;53.70310, $         ; He I
              53.70293, $         ; NIST He I +0.04
              54.38860, $         ; Ne IV +0.02
              55.00318, $         ; Al XI -0.008

              ; the next line is bright but ambiguous from many blends
              ;55.45130, $       ; O IV -0.42 blended with OIII and OIV 55.4076
              
              ;56.27030, $       ; Ne VI +0.31
              56.27980, $       ; Ne VI -0.16
              56.97590, $         ; Ne V +0.11
              ;56.98370, $         ; Ne V -0.28
              57.23360, $         ; Ne V -0.17
              57.40700, $         ; O III -0.09
              ;57.42810, $       ; C III -1.14
              58.09202, $         ; Si XI -0.06
              ;58.09710, $         ; O II 
              ;58.43350, $         ; He I
              58.43339, $         ; NIST He I -0.03
              59.95900, $         ; O III +0.03
              ;60.41212, $         ; Si XI 0.34 blend? Ar III 60.4159, Ne VI 60.43
              ;60.98290, $         ; O IV -0.20 blend with O III 60.9705
              62.49426, $         ; Mg X -0.02
              62.97320, $         ; O V -0.01
              65.73190, $         ; S IV -0.001
              66.13960, $         ; S IV -0.03

              ;68.58170, $       ; N III  blend -0.38
              68.5740, $ ; NIST N III 0.008
              ;70.60614, $         ; Mg IX (weak S II at 70.6144) -0.24
              71.85060, $         ; O II 0.03
              75.0150, $ ; NIST Ne V
              ;75.02210, $       ; S IV -0.40
              ;76.04460, $         ; O V (blend with 76.0227) -0.28
              76.51470, $         ; N IV -0.14
              77.04103, $         ; Ne VIII -0.13
              78.03254, $         ; Ne VIII -0.12
              78.77100, $         ; O IV -0.08
              79.01990, $         ; O IV -0.05
              ;83.3715, $ ; O III
              ;83.37490, $       ; O III (blends with O II and O III) -0.6
              ;83.44670, $         ; O II
              90.4053, $ ; BLEND C II (90.3963+90.4143)/2. 
              ;90.41430, $       ; C II (blend with 90.3963) -0.44
              91.93420, $         ; H I NIST lyman-kappa -0.07
              92.09470, $         ; H I NIST lyman-iota +0.03
              92.31480, $         ; H I NIST lyman-theta -0.09
              92.62290, $         ; H I NIST lyman-eta  -0.16
              ;92.62490, $         ; H I NIST lyman-eta -0.26
              93.07510, $         ; H I NIST lyman-zeta -0.17
              93.33800, $         ; S VI +0.008
              93.78010, $         ; H I NIST lyman-epsilon 4 H I transitions
              ;93.78140, $         ; H I NIST lyman-epsilon -0.21
              94.45240, $         ; S VI (Blend with Ne IV at 94.4405) -0.25
              94.97420, $       ; H I NIST lyman-delta -0.05
              97.25170, $         ; H I NIST lyman-gamma +0.03
              97.70200, $         ; C III -0.03
              98.97990, $         ; N III blend wtih Si II 98.9873 +0.23
              99.15770, $         ; N III -0.03
              ;99.91830, $         ; Ne VI +1.38
              ;99.9497, $        ; NIST O I blended -0.19
              
              101.00850, $      ; C II blend with 101.0373 0.34
              102.57280, $        ; H I NIST lyman-beta -0.19
              103.19138, $        ; O VI +0.02
              103.76154 $         ; O VI  Blend with Fe III at 103.7455-0.36
             ]

  return,reflines
end

function get_eve_ref_lines
  reflines = get_ref_lines()
  gd=lindgen(n_elements(reflines))
  gd=gd[4:*]
  gd=[2,[gd]]
  evereflines = reflines[gd]
  return, evereflines
end


pro find_peaks, reflines, wave, irr, peaksnm, gfitparams, name=name, nostop=nostop, noplot=noplot

  peaksnm = double(reflines)

  searchHalfWidth = 2

  nterms=4 ; gaussfit nterms parameter
  gfitparams = dblarr(n_elements(reflines),nterms)

  if size(nostop,/type) eq 0 then nostop=1
  ;nostop=0 ; default is to stop each iteration
  ;nostop=1 ; default is to not stop

  if keyword_set(noplot) eq 0 then begin
     !P.multi=[0,1,2]
     !p.color=0 & !p.background='ffffff'x & !p.charsize=1.5
  endif
  
  ; calculate the peak locations in the native nm scale wave
  for i=0,n_elements(peaksnm)-1 do begin
     ; locate closest pixel to theoretical peak
     tmp = min(abs(reflines[i] - wave), thisidx )
     ; find local max within 5 bins on either side
     tmp = max(irr[ (thisidx-searchHalfWidth): (thisidx+searchHalfWidth)], maxidx )
     themaxidx = maxidx + thisidx-searchHalfWidth ; relative to index 0
     
     ; TODO: should check that the max is the right peak and not a neighbor
     
     ; fit a quadratic to the peak bins centered at naxidx
     lo = themaxidx - 2
     hi = themaxidx + 2
     ; to avoid "Warning: invert detected a small pivot element"
     ; subtract the center wavelength
     dw = wave[lo:hi] - wave[themaxidx]
     irr_scaled = ts_scale( irr[lo:hi] ) ; 0 to 1

     ; try mpfit to restrict line width to be greater than something
     start = dblarr(nterms)
     start[0] = 1e-7
     start[1] = reflines[i]
     start[2] = 0.04
     start[3] = 1e-7
     pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},4)
     pi[1].limited = 1 & pi[1].limits=reflines[i]+[-.04,.04] ; peak location
     pi[2].limited = 1 & pi[2].limits=[0.031,.042] ; width
     err = sqrt(1./(irr[lo:hi] > (1.d-8)))
     mppar = mpfitfun('mygaussnterms4', wave[lo:hi],irr[lo:hi], err, start, $
                      parinfo=pi,/quiet)
     gfitparams[i,*] = mppar ; could reconstruct just lines later
     
     ; robust polynomial fit
     ; TODO: use weights to weight center the most and edges the least

     r = svdfit(dw, irr_scaled, 3, /double, yfit=yfit)
     
     ; find where the derivative is 0
     ; since dy/dx = r[1] + 2*r[2]*peakwave = 0
     ; then peakwave = -1*r[1] / (2*r[2])
     peaksnm[i] = ( -1.d0 * r[1] / (2.d0 * r[2]) ) + wave[themaxidx]

     if keyword_set(noplot) eq 0 then begin
        gcolor='7700'x
        quadcolor='fe0000'x
        plot,wave[lo-2:hi+2],irr[lo-2:hi+2],xs=1,ys=1,ps=-4, $
             ytit='W/m^2/nm',xtit='Wavelength (nm)', $
             tit=name+' Ref line = '+strtrim(reflines[i],2)+' nm'
        oplot,wave[lo:hi],mygaussnterms4(wave[lo:hi],mppar),ps=-1,co=gcolor,thick=3
        oplot,reflines[i]*[1,1],!y.crange,thick=2
        oplot,peaksnm[i]*[1,1],!y.crange,lines=2,co=quadcolor
        oplot,mppar[1]*[1,1],!y.crange,lines=4,co=gcolor
     
        plot,dw,irr_scaled,xr=[dw[0],dw[-1]],ys=1,ps=-4,$
             ytit='Scaled (arb)',xtit='Relative Wavelength (nm)', $
             tit=name+' Ref='+strtrim(reflines[i],2)+' gausspeak@'+strtrim(mppar[1],2)
        oplot,dw,yfit,co=quadcolor,ps=-1
        oplot,dw,ts_scale(mygaussnterms4(wave[lo:hi],mppar)),ps=-1,co=gcolor
        oplot,reflines[i]*[1,1],!y.crange,lines=1
        oplot,peaksnm[i]*[1,1],!y.crange,lines=1,co=quadcolor
        oplot,mppar[1]*[1,1],!y.crange,lines=1,co=gcolor

        xyouts,-.02,.5,'quadratic shift = '+strtrim(string((peaksnm[i]-reflines[i])/.02,form='(f7.3)'),2)+' bins',co=quadcolor,charsize=1.5
        xyouts,-.02,.6,'gaussian shift = '+strtrim(string((mppar[1]-reflines[i])/.02,form='(f7.3)'),2)+' bins',co=gcolor,charsize=1.5
        xyouts,-.02,.1,'Wavelength peak bin = '+strtrim(string(wave[themaxidx],form='(f7.2)'),2)+' nm'
     endif ; noplot

     if nostop eq 0 then begin
        print,'set nostop=1 to prevent stopping'
        stop
     endif

  endfor

  if keyword_set(noplot) eq 0 then begin
     !p.multi=[0,1,2]
     plot,reflines,(gfitparams[*,1]-reflines)/0.02, ps=-4,xr=[33,105],xs=1,$
          xtit='Reference Line Wavelength (nm)', ytit='fit-ref/0.02 (bins)', $
          yrange=[-1,1],ys=1,tit=name
     plot,reflines,gfitparams[*,2],yr=[.03,.05],ys=1,ps=-4,xr=[33,105],xs=1,$
          xtit='Reference Line Wavelength (nm)', ytit='fit width term', tit=name
     stop
  endif
  
  return
end

;+
; interpolate darr[1,*] to force peakctr-reflines to be zero.
; reflines is a 1d array of reference wavelengths (nm)
; peakctr is the corresponding measurement peak location (nm)
; darr is a 2d array of [3,n_wave] where [0,*] is wavelength sampling
; darr[1,*] is the irradiance and darr[2,*] is an uncertainty
;
; returns a copy of darr with the irradiance shifted to the wavelength grid
;-
function shift_spectrum, reflines, darr, peakctr
  wave = reform(darr[0,*])
  oldirr = reform(darr[1,*])
  new = oldirr ; initialize to same scale

  ; calculate difference of reference lines and measurement peaks

  delta = peakctr-reflines  ; difference
  offset = interpol( delta, reflines, wave )
  ; make first part constant

  short = where(wave le reflines[0])
  offset[short] = 0.d ; use a ramp to delta[0] 
  offset[short[-1]-9:short[-1]] = (dindgen(10)/10.)*delta[0]
  ; make last part contant
  offset[where(wave ge reflines[-1])] = delta[-1]

  newwave = wave + offset

  new = darr ; copy
  new[1,*] = interpol(oldirr, wave, newwave)
  ;stop
  return,new
end

function get_mean_eve, filename

  data = mrdfits(filename,1,/unsign)

  ; fix counts
  data.cps_per_pixel += 0.015
  
  ; need level 1 wavelength scale - copied from megs_wave_defines.pro
  min_megsb_L1_wave = 33.0d0    ;nm
  max_megsb_L1_wave = 107.0d0   ;nm
  num_megsb_L1_spectral_elements  = $
        LONG((max_megsb_L1_wave  - min_megsb_L1_wave) * 50)  ;33.0 - 107.0 nm

  min_megsa1_L1_wave = 3.0d0    ;nm
  max_megsa1_L1_wave = 39.0d0   ;nm (primary 5-19)
  num_megs_L2_spectral_elements = $
        LONG((max_megsb_L1_wave - min_megsa1_L1_wave ) * 50) ; 3 - 107 nm
  megsb_L1_wave  = min_megsb_L1_wave  + (0.02d0)*dindgen(num_megsb_L1_spectral_elements)
  megs_L2_wave = min_megsa1_L1_wave + (0.02d0)*dindgen(num_megs_L2_spectral_elements) ; 3-106.98 inclusive


  ; calculate a spectral mean
  gd=where(data.irradiance[870] > 1e-6) ; need some signal in He continuum
  ; need to filter down to a reference

  sp = median(data[gd].cps_per_pixel,dim=2)
  eve=dblarr(3,n_elements(megsb_l1_wave))
  eve[0,*] = megsb_l1_wave
  eve[1,*] = sp

  return,eve
end

pro adjust_eve

  ; eve is equivalent to rocket 2d array, but only MEGS-B
  evereflines = get_eve_ref_lines()

  ; look at one of the first days, why not the rocket day
  ; compare to EVE L1 (no degradation, no adjustments)
  files = ['MB__L1_2010123_18_007_01.fit.gz', $
           'MB__L1_2014150_08_007_01.fit.gz', $
           'MB__L1_2014360_20_007_01.fit.gz', $
           'MB__L1_2020040_12_007_01.fit.gz', $
           'MB__L1_2021001_15_007_01.fit.gz']

  org=dblarr(n_elements(files),3700) ; MEGS-B
  adj=org
  for i=0,n_elements(files)-1 do begin
     filename = files[i]

     eve = get_mean_eve( filename ) ; unshifted
     org[i,*] = eve[1,*]
     
     find_peaks, evereflines, reform(eve[0,*]), reform(eve[1,*]), peaksnm2010, evegfitparams2010,name='EVE '+file_basename(filename),nostop=1
     print,'found peak locations'
     stop
     ; this is the part needed for production processing
     ; the evegfitparams are the gaussian fits, use the fit center to shift
     eveadj = shift_spectrum( evereflines, eve, evegfitparams2010[*,1] )
     adj[i,*] = eveadj[1,*]
     print,'eveadj calculated'
     stop
  endfor

  ; look for changes within an hour
  name=files[-2]
  hdata = mrdfits(name,1,/unsign) ; next to last one
  
  ; fix counts
  hdata.cps_per_pixel += 0.015

  hep = dblarr(n_elements(hdata), n_elements(evereflines), 4)
  hadj = dblarr(n_elements(hdata), 3700)
  for i=0,n_elements(hdata)-1 do begin
     find_peaks, evereflines, reform(eve[0,*]), reform(hdata[i].cps_per_pixel), hp, hevegfitparams,name='EVE 60sec spectra '+name,/nostop
     hep[i,*,*] = hevegfitparams
     heve=eve

     heve[1,*] = hdata[i].cps_per_pixel
     tmp = shift_spectrum( evereflines, heve, hevegfitparams[*,1] )
     hadj[i,*] = tmp[1,*]
  endfor

  !p.multi=[0,1,1]
  plot,eve[0,*],mean(hdata.cps_per_pixel,dim=2),xr=[33,33.6],xtit='Wavelength (nm)',ytit='L1 cps (arb)'
  for i=0,n_elements(hdata)-1 do oplot,eve[0,*],hdata[i].cps_per_pixel
  stop
  !p.multi=[0,1,2]
  plot,evereflines, (hep[0,*,1]-evereflines)/.02, yr=[-1,1],ys=1,xtit='Reference Line Wavelength (nm)',ytit='Fit-ref/0.02 (bins)',xr=[33,105],xs=1,tit=files[-2]
  for i=0,n_elements(hdata)-1 do oplot,evereflines,(hep[i,*,1]-evereflines)/.02,ps=-4
  eltr = transpose(rebin(evereflines, n_elements(evereflines), n_elements(hep[*,0,0]))) ; expand to same dims as hep (hour EVE fit params)
  delta = (hep[*,*,1]-eltr)/.02
  mdelta = mean(delta, dim=1)
  plot, evereflines, mdelta,ps=-4,ys=1,tit='Mean offset (bins)', xr=[33,105], xs=1
  st = stddev(delta, dim=1)
  for i=0,n_elements(evereflines)-1 do oplot, evereflines[i]*[1,1], mdelta[i]+st[i]*[-1.,1.], co='fe'x
  oplot, !x.crange, [0,0], co='fe00'x, thick=2
  
  
  stop

  !p.multi=[0,1,2]
  ; now plot original spectra
  plot,eve[0,*],mean(hdata.cps_per_pixel,dim=2),xr=[33,33.6],xtit='Wavelength (nm)',ytit='L1 cps (arb)',tit='Native spectra'
  for i=0,n_elements(hdata)-1 do oplot,eve[0,*],hdata[i].cps_per_pixel

  ; now plot adjusted spectra
  plot,eve[0,*],mean(hadj,dim=1),xr=[33,33.6],xtit='Wavelength (nm)',ytit='L1 irrad (arb)',tit='Wavelength adjusted'
  for i=0,n_elements(hdata)-1 do oplot,eve[0,*],hadj[i,*]

  stop
  return
end



pro fix_rocket_wavelengths


  window,0,xs=10,ys=10
  plot,findgen(10)
  wdelete
  
  ; closely examine the wavelengths and identify shifts

  reflines = get_ref_lines()

  d2010=read_dat('36258/rkt_36258_irradiance_v1_1_02nm.dat')
  d2013=read_dat('36290/rkt_36290_irradiance_v1_0_02nm.dat')
  d2018=read_dat('36336/rkt_36336_irradiance_v3_9_02nm.dat')
  d2021=read_dat('36353/code/rkt_36353_irradiance_v1_0_02nm.dat')
  d2021[1,*] = d2021[1,*] > (d2018[1,*]*1.01)
  
  !p.multi=[0,1]
  !p.color=0 & !p.background='ffffff'x & !p.charsize=1.5
  plot,d2010[0,*],d2010[1,*],xr=[33,106],xs=1,yr=[1e-7,1e-3],ys=1,/ylog,$
       xtit='Wavelength (nm)',ytit='W/m^2/nm',tit='2010,2013,2018,2021 rocket spectra'
  oplot,d2013[0,*],d2013[1,*],co='fe'x
  oplot,d2018[0,*],d2018[1,*],co='fe0000'x
  oplot,d2021[0,*],d2021[1,*],co='fe00'x
  for i=0,n_elements(reflines)-1 do oplot,reflines[i]*[1,1],[1e-7,1e-3],co='aaaa'x
  stop

  plot,d2010[0,*],d2010[1,*],xr=[33,106],xs=1,yr=[1e-7,1e-3],ys=1,/ylog,$
       xtit='Wavelength (nm)',ytit='W/m^2/nm',tit='2010,2013,2018,2021 rocket spectra'
  oplot,d2013[0,*],d2013[1,*],co='fe'x
  oplot,d2018[0,*],d2018[1,*],co='fe0000'x
  oplot,d2021[0,*],d2021[1,*],co='fe00'x
  stop
  
  ; use linear interpolation to move peaks to match reference lines
  ; quadratic fit to top 4 bins, peak is where derivative=0

  find_peaks, reflines, reform(d2010[0,*]), reform(d2010[1,*]), peaksnm2010, gfitparams2010,name='2010 36.258',/nostop,/noplot
  find_peaks, reflines, reform(d2013[0,*]), reform(d2013[1,*]), peaksnm2013, gfitparams2013,name='2013 36.290',/nostop,/noplot
  find_peaks, reflines, reform(d2018[0,*]), reform(d2018[1,*]), peaksnm2018, gfitparams2018,name='2018 36.336',/nostop,/noplot
  find_peaks, reflines, reform(d2021[0,*]), reform(d2021[1,*]), peaksnm2021, gfitparams2021,name='2021 36.353',/nostop,/noplot

  
  ; perform a linear shift between peaknm
  newd2010 = shift_spectrum( reflines, d2010, gfitparams2010[*,1] )
  newd2013 = shift_spectrum( reflines, d2013, gfitparams2013[*,1] )
  newd2018 = shift_spectrum( reflines, d2018, gfitparams2018[*,1] )
  newd2021 = shift_spectrum( reflines, d2021, gfitparams2021[*,1] )

  ; store new rocket spectra
  form='(f7.3,e13.6,f7.3)'
  coltext = 'wavelength (nm), irradiance (W/m^2/nm), relative uncertainty (fraction)'
  write_dat, file='rkt_36258_irradiance_v1_1_02nm_shifted.dat', newd2010,form=form,col=coltext
  write_dat, file='rkt_36290_irradiance_v1_0_02nm_shifted.dat', newd2013,form=form,col=coltext
  write_dat, file='rkt_36336_irradiance_v3_9_02nm_shifted.dat', newd2018,form=form,col=coltext
  write_dat, file='rkt_36353_irradiance_v1_0_02nm_shifted.dat', newd2021,form=form,col=coltext
  print,'wrote new rocket spectrum dat files'
  
  !p.multi=[0,1,1]
  ; full range
  plot,newd2010[0,*],newd2010[1,*],xr=[33,106],xs=1,yr=[1e-7,1e-3],ys=1,/ylog,$
       xtit='Wavelength (nm)',ytit='W/m^2/nm',tit='2010,2013,2018,2021 tuned rocket spectra'
  oplot,newd2013[0,*],newd2013[1,*],co='fe'x
  oplot,newd2018[0,*],newd2018[1,*],co='fe0000'x
  oplot,newd2021[0,*],newd2021[1,*],co='fe00'x
  stop
  
  !p.multi=[0,1,2]
  ; loop over 2 nm chunks to view line centers
  minsp = d2010[1,*]<d2013[1,*]<d2018[1,*]<d2021[1,*]
  maxsp = d2010[1,*]>d2013[1,*]>d2018[1,*]>d2021[1,*]
  for lownm = 33,103,2 do begin
     gd=where(d2010[0,*] ge lownm and d2010[0,*] lt lownm+2.1)
     yr = [min(minsp[gd]),max(maxsp[gd])]
     plot,d2010[0,*],d2010[1,*], xr=[lownm,lownm+2.1], xs=1, tit='Rocket native wavelengths',yr=yr,ys=1,ps=10
     oplot,d2013[0,*],d2013[1,*], ps=10, co='fe'x
     oplot,d2018[0,*],d2018[1,*], ps=10, co='fe0000'x
     oplot,d2021[0,*],d2021[1,*], ps=10, co='fe00'x
     for i=0,n_elements(reflines)-1 do oplot,[1,1]*reflines[i],!y.crange,lines=1

     plot,newd2010[0,*],newd2010[1,*], xr=[lownm,lownm+2.1], xs=1,tit='Rocket tuned wavelength', xtit='Wavelength (nm)',yr=yr,ys=1, ps=10
     oplot,newd2013[0,*],newd2013[1,*], ps=10, co='fe'x
     oplot,newd2018[0,*],newd2018[1,*], ps=10, co='fe0000'x
     oplot,newd2021[0,*],newd2021[1,*], ps=10, co='fe00'x
     for i=0,n_elements(reflines)-1 do oplot,[1,1]*reflines[i],!y.crange,lines=1
     
     stop
  endfor
  
  stop
  return
end

  
