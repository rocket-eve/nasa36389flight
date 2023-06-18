function continuum_fit, wave, irradiance
  ; function is a line in log-log space
  coeffs = poly_fit( alog(wave), alog(irradiance), 1)
  return, coeffs
end

;+
; This is the function for use with mpfit to fit
; the Planck function from Machado et al 2018 to return
; temperature and b.
; Here P has 2 elements [ b, T ]
; b is dimensionless, T is Kelvin
;-
function planck_function, wave, P
  sunradius=6.957d8 ; meters (https://arxiv.org/pdf/1510.07674.pdf)
  area = !dpi * sunradius * sunradius ; m^2 ~1.52e18
  au = 1.495978707d11 ; meters (10.1007/s10569-009-9203-8)
  a = area/(100.d*au*au)        ; no units ~6.79e-7
  ; 100ergs/cm^2/sec/Angstrom=1W/m^2/nm
  h = 6.62607015d-34 ; Plancks constant, Watts*sec*sec
  c = 2.99792458d8 ; meters/sec
  kb = 1.380649d-23 ; Boltzmans constant Joule/Kelvin or Watts*sec/Kelvin (SI base unit)
  twohcc = 2.d * h * c * c ; ~1.19e-16
  hcoverkb = h * c / kb ; ~0.0144

  lambda = wave * 1.d-9 ; meters
  
  ; b=P[0] must be greater than 0.1 and less than 10000 for H LyC
  ; T=P[1] must be greater than 3000 and less than 1e5 for H LyC
  
  ;f = (a / P[0]) * (1.d-9) * (twohcc / (lambda^5)) / $
  ;    (exp( hcoverkb / (lambda * P[1])) - 1.d0)

  n1 = a*(1.d-9)*twohcc/(lambda^5) ; ~8.09e13 to 5.26e3 at 1nm and 110 nm
  d1 = hcoverkb/lambda ; ~1.44e7 to 1.32e5 at 1nm and 110 nm
  f = (n1/P[0]) / (exp( d1 / P[1] ) - 1.d0)
  
  return, f
end


pro fit_h_cont, wave, irradiance, hydrogencoeffs=hydrogencoeffs, helium1coeffs=helium1coeffs, carbon1coeffs=carbon1coeffs

  no_plots = 0 ; set to 0 to show plots, 1 to suppress plots
  
  if size(wave,/type) eq 0 or size(irradiance,/type) eq 0 then begin
     @config36389
     !p.charsize=1.5
     !x.margin=[11,3]
     ;rocketspectrum = read_dat('rkt_36389_irradiance_v1_0.dat')
     ;wave = reform(rocketspectrum[0,*])
     ;irradiance = reform(rocketspectrum[1,*])
     ;irr_uncertainty = reform(rocketspectrum[2,*])
     
     workingdir=file_dirname(routine_filepath()) ; in /code/
     datadir = file_dirname(workingdir)+'/data/'

     ;restore,datadir+'/mb_36389_irr_at_1au.sav' ; wave, irradiance, relEerr

     ; try another rocket
     ;restore,datadir+'/mb_36353_irr_at_1au.sav' ; wave, irradiance, relEerr
     restore,'../../36336/rkt_36336_irradiance.sav' ; rkt
     wave = reform(rkt[0,*])
     irradiance = reform(rkt[1,*])
     relEerr = reform(rkt[2,*])
     
     
     irr_uncertainty = relEerr
     ;stop
  endif
  
  h_gd=where( $
       (wave gt 79.2 and wave lt 83.13) or $ ; O II/O III
       (wave gt 83.7 and wave lt 85.7) or $ ; cut a hole at ~85.804 C II/Mg VII/S V
       (wave gt 86.0 and wave lt 89.4) or $ ; cut a hole at ~89.517 Ne VII
       (wave gt 89.6 and wave lt 90.2) or $ ; cut a hole for C II ~90.414
       (wave gt 90.6 and wave lt 91.12)) ; -67.730147, 12.887809

  ; create a fit
  hydrogencoeffs = continuum_fit( wave[h_gd], irradiance[h_gd] )
  print,'INFO: H continuum coefficients:',strtrim(hydrogencoeffs,2)

  if ~no_plots then begin
     !p.multi=[0,1,2]
  
     plot,wave,irradiance,/ylog,xr=[78,92],xs=1,yr=[5e-7,1e-3],ys=1, $
          ytit='Irradiance', xtit='Wavelength (nm)'
     oplot,wave[h_gd],irradiance[h_gd],co='fe'x,ps=4
     oplot, wave, exp(poly(alog(wave), hydrogencoeffs)),co='fe'x
  
     ; the coefficients are related to the plasma temperature
     ; somehow...
     ; maybe there is a way to compare with some synthetic spectra
     ; to infer what the temperature is from those  
  
     hconteval = exp(poly(alog(wave), hydrogencoeffs))
     bad=where(wave gt 91.1267) ; ~ 
     hconteval[bad]=0.          ; no hydrogen continuum longer than 
  
;     stop
     ;try to fit using the plank function to get b and T
     start_params = [10.4d0, 9300.d0]
     parinfo = replicate({value:0.d, fixed:0, limited:[1,1], limits:[0.d,0.d]},n_elements(start_params))
     parinfo[0].limits=[0.5d, 2000.d]
     parinfo[1].limits=[3000.d, 20000.d]
     ;err = sqrt(irradiance[h_gd]) ; too small
     ;err = irradiance[h_gd]*.01 ; 20% 2.3,113.3, 10% 1.18, 56.6, 1% 0.118,5.66
     err = irr_uncertainty[h_gd]*irradiance[h_gd] ; use total uncertainty
     params = mpfitfun('planck_function', wave[h_gd], irradiance[h_gd], err, start_params,parinfo=parinfo, perror=perror, ftol=1d-30, covar=covar)
     oplot,wave,planck_function(wave,params),co='fe00'x

     plot,wave[h_gd], hconteval[h_gd] - irradiance[h_gd],ys=1,ps=-4,ytit='Diff (W/m^2/nm)',$
          xrange=!x.crange, xs=1, $
          tit='Fit minus measurement',xtit='Wavelength (nm)'
     oplot,wave[h_gd], hconteval[h_gd]-irradiance[h_gd],co='fe'x,ps=-4
     oplot,wave[h_gd], planck_function(wave[h_gd], params)-irradiance[h_gd],co='fe00'x,ps=-4
     print,'H LyC Temperature = '+strtrim(params[1],2)+' +/- '+strtrim(perror[1],2)+' logT='+strtrim(alog10(params[1]),2)
     print,'H LyC b (NLTE) = '+strtrim(params[0],2)+' +/- '+strtrim(perror[0],2)
     
     stop
     
  endif

  ; try to fit the helium 1 continuum
  he_gd=where($
       (wave gt 47.1 and wave lt 47.9) or $
       (wave gt 49.25 and wave lt 49.7) or $
       (wave gt 50.15 and wave lt 50.4) $
            )

  helium1coeffs = continuum_fit( wave[he_gd], irradiance[he_gd] )
  print,'INFO: He I continuum coefficients:',strtrim(helium1coeffs,2)

  if ~no_plots then begin
     plot,wave,irradiance,/ylog,xr=[44,52],xs=1,yr=[5e-7,1e-3],ys=1, $
          ytit='Irradiance', xtit='Wavelength (nm)'
     oplot,wave[he_gd],irradiance[he_gd],co='fe'x,ps=4
  
     heconteval = exp(poly(alog(wave), helium1coeffs))
     bad = where(wave gt 50.42) ; estimated by 24.5874eV (1239.838/eV
     heconteval[bad] = 0.       ; no hydrogen continuum longer than 

     oplot,wave,heconteval,co='fe'x ; red line

     ;try to fit using the plank function to get b and T
     start_params = [10.4d0, 20300.d0]
     parinfo = replicate({value:0.d, fixed:0, limited:[1,1], limits:[0.d,0.d]},n_elements(start_params))
     parinfo[0].limits=[0.001d, 2000.d]
     parinfo[1].limits=[3000.d, 200000.d]

     err = irr_uncertainty[he_gd]*irradiance[he_gd] ; use total uncertainty
     params = mpfitfun('planck_function', wave[he_gd], irradiance[he_gd], err, start_params,parinfo=parinfo, perror=perror, ftol=1d-30, covar=covar)
     oplot,wave,planck_function(wave,params),co='fe00'x

     
     plot,wave[he_gd], heconteval[he_gd] - irradiance[he_gd],ys=1,ps=-4,$
          ytit='Diff (W/m^2/nm)', xr=!x.crange, xs=1,$
          tit='Fit minus measurement',xtit='Wavelength (nm)'
     
     oplot,wave[he_gd], heconteval[he_gd] - irradiance[he_gd],ps=-4,co='fe'x
     oplot,wave[he_gd], planck_function(wave[he_gd], params)-irradiance[he_gd],co='fe00'x,ps=-4
     print,'He LyC Temperature = '+strtrim(params[1],2)+' +/- '+strtrim(perror[1],2)+' logT='+strtrim(alog10(params[1]),2)
     print,'He LyC b (NLTE) = '+strtrim(params[0],2)+' +/- '+strtrim(perror[0],2)

     stop
  endif

  ; carbon continuum
  c_gd = where( $
         (wave gt 96.0 and wave lt 97.0) or $ ; very dim
         (wave gt 100.2 and wave lt 100.5) or $
         (wave gt 104.25 and wave lt 104.5) or $
;         (wave gt 104.65 and wave lt 104.71) or $
;         (wave gt 105.09 and wave lt 105.13) or $
         (wave gt 104.65 and wave lt 105.13) or $
         (wave gt 105.5 and wave lt 105.79) ) ; only available in sav file
  carbon1coeffs = continuum_fit( wave[c_gd], irradiance[c_gd] )
  print,'INFO: C I continuum coefficients:',strtrim(carbon1coeffs,2)

  if ~no_plots then begin
  
     ;xr=[92,105.8]
     xr=[long(wave[c_gd[0]]), ceil(wave[c_gd[-1]])]
     ;yr=[5e-7,1e-3]
     yr=[1e-6,2e-5]
     plot,wave,irradiance,/ylog,xr=xr,xs=1,yr=yr,ys=1, $
          ytit='Irradiance', xtit='Wavelength (nm)'
     oplot,wave[c_gd],irradiance[c_gd],co='fe'x,ps=4
  
     cconteval = exp(poly(alog(wave), carbon1coeffs))
     ; carbon 1 continuum stops at 110.1 nm
     bad = where(wave gt 110.1, n_bad) ; estimated by 11.2603eV (1239.838/eV
     if n_bad gt 0 then cconteval[bad] = 0.       ; no c continuum longer than that

     oplot,wave,cconteval,co='fe'x ; red line


     ;try to fit using the plank function to get b and T
     start_params = [2.6d0, 7000.d0]
     parinfo = replicate({value:0.d, fixed:0, limited:[1,1], limits:[0.d,0.d]},n_elements(start_params))
     parinfo[0].limits=[0.5d, 2000.d]
     parinfo[1].limits=[3000.d, 20000.d]

     err = irr_uncertainty[c_gd]*irradiance[c_gd] ; use total uncertainty
     params = mpfitfun('planck_function', wave[c_gd], irradiance[c_gd], err, start_params,parinfo=parinfo, perror=perror, ftol=1d-30, covar=covar)
     oplot,wave,planck_function(wave,params),co='fe00'x
     
     plot,wave[c_gd], cconteval[c_gd] - irradiance[c_gd],ys=1,ps=-4,$
          ytit='Diff (W/m^2/nm)', xr=xr,xs=1, $
          tit='Fit minus measurement',xtit='Wavelength (nm)'
     oplot,wave[c_gd], cconteval[c_gd] - irradiance[c_gd],ps=-4,co='fe'x
     oplot,wave[c_gd], planck_function(wave[c_gd], params)-irradiance[c_gd],co='fe00'x,ps=-4
     print,'C LyC Temperature = '+strtrim(params[1],2)+' +/- '+strtrim(perror[1],2)+' logT='+strtrim(alog10(params[1]),2)
     print,'C LyC b (NLTE) = '+strtrim(params[0],2)+' +/- '+strtrim(perror[0],2)

     stop
  endif
  stop
  
  return
end
