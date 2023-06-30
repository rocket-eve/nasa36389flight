function fitfunc, X, P
  f = p[0] + p[1]*exp( p[2]*x )
  return,f
end

;+
; This function calculates the curve fit of y vs x with err representing
; uncertainties in y (sqrt(y)).
;
; :Params:
;    x: in, required, type=fltarr
;      X-values (time) for each point to fit
;    y: in, required, type=fltarr
;      Y-values (cps) for each point to fit
;    err: in required, type=fltarr
;      error in each y value, use sqrt(y) if not known
;    perror: out, optional, type=floatarr
;      The uncertainty in the fit coefficients.
;    
;-
function get_fit_params, x, y, err, perror, widx=widx, wave=wave

  maxy = max(y,min=miny)
  lastquartilemedian = median(y[n_elements(y)/4 : *])

  if maxy lt 1 then begin
     perror=[0.,0.,0.]
     cpsfit=[0.,0.,0.]
     return,cpsfit
  endif
  
  start = [maxy, -14.*maxy, -1e-3]
  ;  y = p0 + p1*exp(p2*t)
  parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
                       limits:[0.D,0]}, n_elements(start))
  parinfo[0].limited[*]=1
  parinfo[0].limits = lastquartilemedian * [0.9, 1.3] ; limit max absorption to 30%*max
  parinfo[0].limits[1] = (maxy-miny)*0.5 + maxy
  parinfo[1].limited[1] = 1
  parinfo[1].limits[1] = -0.5 * maxy
  parinfo[2].limited[1] = 1
  parinfo[2].limits[1] = -1e-4
     
  cpsfit = mpfitfun('fitfunc', x, y, err, start, $
                    parinfo=parinfo, $
                    perror=perror, covar=covar, $
                    status=status, $
                    /quiet)
  ; when the fit coefficients reach the limit, the uncertainty that
  ; is returned is 0. Replace zero with a more reasonable estimate.
  ; Here it is set to 3*(RMS difference from the fit)
  perror[0] = perror[0] > 3.*(sqrt(mean((fitfunc(x,cpsfit) - y)^2)))

  ;plot,x,y,ps=-4,xs=1,ys=1,xtit='Time since launch (seconds)', $
  ;     ytit='Counts @ '+strtrim(wave,2)+'nm',$
  ;     tit='Index=#'+strtrim(widx,2)+' extrap='+strtrim(cpsfit[0],2)+'+/-'+strtrim(perror[0],2)
  ;oplot,x,fitfunc(x,cpsfit),co='fe'x
  ;stop

  
  return,cpsfit
end


pro study_lambda, spmodel, slit1=slit1, slit2=slit2, megsb=megsb

  pm_in = !p.multi
  pc_in = !p.charsize
  
  ; use gaussfit to get a gaussian for 36.0 and 36.8 nm
  n_sp = n_elements(spmodel.meas)

  if keyword_set(slit2) then begin
     wavelength = double(spmodel.wave2)
     ; MEGS-A slit 2 lines
     targetlinew=[ 14.83770d0, 15.21510d0, 15.41620d0, $ ;15.43030d0, $
                   15.77290d0, 17.10730d0, 17.45310d0, $
                   17.72400d0, 18.04010d0, 18.21670d0, $
                   18.45370d0, 18.68870d0, 18.82990d0, $
                   19.23940d0, 19.51190d0, 20.20440d0, 20.38260d0, $
                   21.13170d0, 21.71010d0, 22.47260d0, 22.70000d0, $
                   23.73310d0, 23.85700d0, 25.19520d0, 25.63170d0, $
                   25.83740d0, 26.10560d0, 26.29760d0, 26.47880d0, $
                   27.05200d0, 27.21270d0, 28.41630d0, 28.84340d0, $
                   28.84340d0, 29.28090d0, 29.61170d0, 29.960990d0,  30.37800d0, $
                   31.62180d0, 31.98390d0, 32.84480d0, 33.54090d0, $
                   34.53120d0, 34.74020d0, 34.98730d0, 35.60370d0, $
                   35.86220d0, 35.96440d0,  $
                   36.07580d0, 36.44670d0, 36.80710d0, 37.38030d0 $
                 ]
  endif
  if keyword_set(slit1) then begin
     wavelength = double(spmodel.wave2)
     ; MEGS-A slit 1 lines
     targetlinew=[ $
                 6.28710d0, $ ; Fe XVI blend with 6.3152 Mg X and/or 6.32953 Mg X
                 6.59651d0, $ ; Fe XVII - 6.69099 is too dim to use

                 6.96820d0, $   ; Fe XV good FWHM, gooe location
                 7.20290d0, $ ; S VII is blended with 7.23127 Mg IX - large FWHM
                 ;7.30920d0, $ ;Ar XVI blend with 7.34720d0 Fe XV
                 ;7.59372d0, $ ; O VIII blend with 7.61968 Fe XXIII - good FWFM
                 7.61968d0, $ ; Fe XXIII blend with 7.59372 O VIII - good FWFM

                 ;8.04890d0, $   ; Fe XX noisy, blend with 8.04539 Ne IX
                 8.27597d0, $ ; Ne IX some noise, narrow
                 8.69990d0, $; Fe XIX blend with 86.9144 Si VII

                 9.19288d0, $   ; Fe XVII blend with 91.9888 Fe XVIII
                 ;9.37810d0, $
                 9.39320d0, $   ; Fe XVIII possible blend with 93.7811 Fe XX and 94.012 Fe X
                 9.61210d0, $ ; Fe X possible blend with 96.0662 Fe XVII and 96.4400 Si V
                 ;9.82601d0, $
                 9.83552d0, $ ; Fe XX blend of 9.83234 Fe XVII and 9.83931 Fe XXI, looks clean but some fits do not converge
                 10.08087d0, $  ; Fe XVII
                 10.37596d0, $  ; Fe XXI blend with Fe XXI 10.38331 Ni XXIV 10.36755 - good width
                 10.53859d0, $ ; Fe XVII
                 ;10.85080d0, $ ; 
                 10.87531d0, $ ;  Mn XVIII blend with Fe XVIII 10.86501
                 11.16947d0, $ ; Fe XIX, blend with Fe XVIII 11.15790
                 11.41849d0, $ ; Fe XVIII blend with 11.44090 Fe XXII
                 11.69099d0, $ ; Fe XVIII blend with Fe XXII 11.71543
                 ; 12.29720d0, 12.38312d0, $ too blended, need to fit both together
                 12.78367d0, $ ; Fe XIX blend with Fe XIX 12.84062 and Ne VII 12.76660
                 13.12400d0, $ ; Fe VIII blend with Si XII 13.15515
                 13.83680d0, $ ; Mn XX blend with Ne VI 13.8387 and Fe XX 13.84789m - noisy
                 14.83770d0, 15.21510d0, 15.41620d0, $ ;15.43030d0, $
                 15.77290d0, 17.10730d0, 17.45310d0, $
                 17.72400d0, 18.04010d0, 18.21670d0, $
                 18.45370d0, 18.68870d0, 18.82990d0, $
                 19.23940d0, 19.51190d0, 20.20440d0, 20.38260d0, $
                 21.13170d0, 21.71010d0 $
                 ]
  endif

  ; MEGS-B
  if size(wavelength,/type) eq 0 then begin
     wavelength = double(spmodel.wavemb)
     ; MEGS-A slit 2 lines
     targetlinew=[ 33.54090d0, $ ; Fe XVI
                   ;34.53120d0, $
                   ; 34.74020d0, $ possible, Si X blend with Fe XVII 34.7816
                   ; 34.98730d0, ; Si IX, blend with Fe XXII 34.9303
                   35.60370d0, $ ; Si X
                   ;35.84740d0, 35.96440d0,  $ cannot see in MEGS-B
                   36.07580d0, $ ; Fe XVI blend with Mn XV 36.0987
                   36.44670d0, $ ; Fe XII
                   36.80710d0, $ ; Mg IX 
                   ;37.38030d0, $ ; some bad stuff at apogee, lots of O III lines blended
                   39.24330d0, $ ; Al IX
                   40.19410d0, $ ; Ne VI blende with Ne VI at 40.1146
                   40.33070d0, $ ; Ne VI blended with Ne VI 40.32600 
                   41.14438d0, $ ; Fe XVIII blend with Na VIII 41.11660
                   41.72580d0, $ ; Fe XVIII blend with S XIV 41.7660
                   44.57000d0, $ ; S XIV
                   41.97140d0, $ ; C IV
                   43.04540d0, $ ; Mg VIII blend with Mg VII 43.11940
                   43.49800d0, $ ; O III blend with Mg VII 43.49230
                   43.67330d0, $ ; Mg VIII possible blend with Fe XIX 43.67545
                   43.91772d0, $ ; Mg IX 
                   44.39737d0, $ ; Mg IX blend with Fe XIV 44.42210
                   44.57011d0, $ ; S XIV possible blend with Fe XIX 44.57638
                   44.73570d0, $ ; Fe XIV
                   45.96270d0, $ ; C III blend with CIII 45.9633, 45.9516
                   46.52210d0, $ ; Ne VII
                   46.74290d0, $ ; Fe XIV blend with wing of 46.5221
                   46.98250d0, $ ; Ne IV blend with Ne IV 46.9875
                   48.14890d0, $ ; Fe XV (He cont) blend Ne V 48.1374
                   48.29940d0, $ ; Ne V (He cont)
                   48.52500d0, $ ; S III (He cont)
                   48.96290d0, $ ; Ne III blend with He cont and Ne III 48.949,49.029
                   49.94060d0, $ ; Si XII (with He I continuum)
                   50.81780d0, $ ; O III blend with O III 50.768 (He cont ends)
                   ; 51.04? unidentified
                   ; 51.22? unidentified
                   51.56170d0, $ ; He I
                   52.06650d0, $ ; Si XII (Fe XVII at 51.97001) blend with dim Ne IV 52.1739, 52.1824
                   52.57940d0, $ ; O III in wing of 52.5794
                   53.70300d0, $ ; He I (O II 53.7832)
                   ;54.13375d0, $ ; Ne IV blend with Fe XX
                   54.20760d0, $ ; Ne IV
                   54.38860d0, $ ; Ne IV significant unknown blend
                   55.00310d0, $ ; Al XI
                   ; 
                   55.45140d0, $ ; O IV blend of O II and O III 
                   55.77660d0, $ ; Ca X blend with Ne VI 55.86
                   56.28050d0, $ ; Ne VI, blue wing has Ne VII 56.173
                   ;56.787, $ ; Fe XX blended heavily (not usually detectable)
                   56.84240d0, $ ; Ne V
                   56.98280d0, $ ; Ne V blend with Ne V 56.9756
                   57.23350d0, $ ; Ne V blend with 57.2105
                   57.42810d0, $ ; C III blend with O III 57.406
                   58.09700d0, $ ; O II blend with Si XI 58.09202 and O II 58.041
                   58.43340d0, $ ; He I blend red wing with C III 58.5428
                   ;59.22350d0, $ ; Fe XIX (not usually detectable), blend with Ar IV 59.3341
                   59.50220d0, $ ; C II weak
                   ;59.69260d0, $ ; Ca VII weaker blend with O III
                   59.78140d0, $ ; O III blue wing blend with Ca VII
                   59.95900d0, $ ; O III blend with weaker O III 59.7814
                   60.41212d0, $ ; Si XI
                   ;60.83970d0, $ ; O IV blend into blue wing of brighter 60.98 Mg IX/O IV blend
                   60.97930d0, $ ; Mg X blend with O IV at 60.983 & maybe O III 61.0745
                   ;61.63040d0, $ ; O II blend with O IV 61.6952, O II 61.7063
                   ; 62.12? unidentified
                   62.49410d0, $ ; Mg X blend with O IV 62.5129, 62.5853 and 62.4619
                   62.97320d0, $ ; O V
                   ; 63.76? and 63.92? unidentified
                   ; 64.22? unidentified
                   ;64.41530d0, $ ; O II blend with O II 64.4162, N II 64.4837, 64.5178
                   ; 64.94? unidentified
                   65.73190d0, $ ; Si IV blend S V 65.8253
                   66.14550d0, $ ; Si IV
                   66.31260d0, $ ; S V
                   67.13860d0, $ ; N II blend with N II
                    ; 67.66? unidentified
                   ;68.06760d0, $ ; S III blend with brighter 68.1488/68.1578
                   68.15780d0, $ ; S III
                   68.58170d0, $ ; N III blend with N III 68.515,68.6336
                   ;68.73450d0, $ ; C II blend with brigher N III above
                   69.05210d0, $ ; C III blenw with N III 69.1193
                   69.46210d0, $ ; Ni XX blend with Na IN 69.4147
                   ;69.66230d0, $ ; S V weak
                   70.01480d0, $ ; S III weak
                   ;70.23370d0, $ ; O III blend with S III 70.2778
                   70.29000d0, $ ; O III blend with 70.3854
                   70.38540d0, $ ; O III blend on blue side with O III
                   70.60365d0, $ ; Mg IX blend with S VI 70.6471
                   71.85670d0, $ ; O II blend with O II
                   ;72.15580d0, $ ; Fe XX no fit (H cont starting)
                   75.02210d0, $ ; S IV (H cont)
                   ;75.86770d0, $ ; O V  (H cont)
                   76.04460d0, $ ; O V blend with O V 76.0227
                   76.51520d0, $ ; N IV
                                ;76.56850d0, $ ; S II
                   77.04103d0, $ ; Ne VII blend N III 77.1545
                   78.03850d0, $ ; Ne VIII
                   78.64680d0, $ ; S V (Blend Fe XXI & H cont)
                   78.77100d0, $ ; O IV
                   79.02010d0, $ ; O IV
                   83.37490d0, $ ; O III blend with O II 83.34465
                   83.52890d0, $ ; O III blend with O III 83.5092, 83.6595 (H cont)
                   ;84.55715d0, $ ; Fe XXII not present
                   90.41410d0, $ ; C II (end of H cont)
                   ; 91.86? 91.96? unidentified
                   92.0947d0, $ ; H I (NIST)
                   ;92.13650d0, $ ; O IV (Same as NIST)
                   92.32250d0, $ ; N IV (H I 92.3151
                   92.6249d0, $ ; H I (NIST)
                   93.07510d0, $ ; H I (NIST)
                   93.33380d0, $ ; S VI
                   93.78140d0, $ ; H I (NIST)
                   94.45230d0, $ ; S VI
                   94.97430d0, $ ; H I
                   ; 95.88? unidentified
                   97.25170d0, $ ; H I (NIST)
                   97.70200d0, $ ; C III
                   98.97990d0, $ ; N III
                   ;99.0190d0, $ ; O I Sumer Atlas
                   ;99.15110d0, $ ; N III
                   99.15770d0, $ ; N III
                   99.47900d0, $ ; Si III
                   99.73870d0, $ ; Si III
                   ;99.92330d0, $ ; Ne VI
                   99.9497d0, $ ; O I (NIST)
                   100.6091d0, $ ; S II
                   ;100.9857d0, $ ; C II
                   101.00710d0, $ ; C II
                   ;101.0369d0, $ ; C II
                   101.2492d0, $ ; S III very weak
                   101.5775d0, $ ; S III
                   ; 101.76? unidentified
                   102.1321d0, $ ; S III
                   102.5722d0, $ ; H I
                   ; 102.76? unidentified
                   103.1912d0, $ ; O VI
                   ;103.6337d0, $ ; C II questionable, on the wing of 103.76
                   103.7613d0 $ ; O VI
                   ;105.8678d0, $ ; Ne IV (NIST)
                   ;105.953d0 $  ; S I (NIST)
                 ]
     ;stop
  endif
  
  !p.multi=[0,4,3]
  !p.charsize=2
  to_fwhm = 2.d*sqrt(2.d*alog(2.d)) ; refer to gaussfit help 
  dlambda = wavelength[1]-wavelength[0] ; wavelength spacing between samples
  wavectr = wavelength ;+ (dlambda*0.5d0) ; wavelength at center of bin
  ; shifting to the center is NOT necessary for the rocket spectra
  ; the wavelength scale is fit using these as if they are the centers

  ; comparison to the flight EVE requires accounting for the offset in SDO, not the rocket
  
  for lineidx=0,n_elements(targetlinew)-1 do begin
      ; loop from apogee down so yrange shows all spectra

     newplot=1
     widths = fltarr(n_sp) ; track widths over time/alt
     centers = widths ; fitted line center
     offsets = widths ; background offset from fit
     times = fltarr(n_sp)
     
     for i=n_sp-1,0,-1 do begin
        junk = min(abs(wavectr-targetlinew[lineidx]),idx)

        gd = idx - 4L + lindgen(9) ; MEGS-A
        if keyword_set(megsb) then gd = idx - 2L + lindgen(5) ; MEGS-B

        case 1 of
           keyword_set(slit2):  thissp = spmodel.meas[i].sp2
           keyword_set(slit1):  thissp = spmodel.meas[i].sp1
           else: thissp = spmodel.meas[i].mbsp
        endcase

        maxy = max(thissp[gd])

        estimates = [maxy, $                 ; P0
                     targetlinew[lineidx], $ ; P1
                     0.1/to_fwhm, $          ; P2
                     maxy*.01]               ; P3


        fity = gaussfit(wavectr[gd], thissp[gd], params,nterm=4,estimates=estimates, measure_errors=sqrt(thissp[gd])>1., yerror=yerror)
        ;stop ; check convergence
        widths[i] = params[2] * to_fwhm
        centers[i] = params[1]
        offsets[i] = params[3]
        
        print,'T='+strtrim(round(spmodel.meas[i].time),2)+' params: ',strjoin(strtrim(params,2),' ')+' fwhm='+strtrim(widths[i],2)+' err='+strtrim(yerror,2)
        if newplot eq 1 then begin
           title=strtrim(targetlinew[lineidx],2)+' nm'
           case 1 of
              keyword_set(slit1) : title='MA1 '+title
              keyword_set(slit2) : title='MA2 '+title
              else: title='MB '+title
           endcase
           plot,wavectr,thissp,ps=-4,xs=1, /ynozero, $ ;ys=1, $
                tit=title, $
                xtit='Wavelength (nm)', ytit='cps',$
                xrange=[wavectr[gd[0]]-.05, wavectr[gd[-1]]+0.05]
        endif else begin
           oplot,wavectr,thissp,ps=-4
        endelse
        oplot,wavectr[gd],fity,co='fe'x,ps=-1 ; just the points in the fit

        newplot=0
     endfor

     ; write the fwhm and width on the spectrum
     xyouts,/data,params[1],!y.crange[0]+.2*(!y.crange[1]-!y.crange[0]),'fwhm='+strtrim(string(median(widths),form='(f5.3)'),2)+'!Cloc='+strtrim(string(median(centers),form='(f8.3)'),2),charsize=1.25,align=0.5

     ; plot width
     xr=[100,300]
     delta = max(abs(widths-0.1)) < 0.2 ; 0.1 is optical resolution spec
     wyr = 0.1 + delta*[-1.,1.] > 0.
     plot,spmodel.meas.time, widths, xs=1,xr=xr,yr=wyr,ys=1,ps=-4, $
          xtit='Time (sec)',ytit='FWHM (nm)',$
          tit=strtrim(targetlinew[lineidx],2)+' nm'
     narrow=where(widths le 0.1,n_narrow,comp=wide,ncomp=n_wide)
     if n_narrow gt 1 then oplot,spmodel.meas[narrow].time,widths[narrow],ps=6,co='aa00'x
     if n_wide gt 1 then oplot,spmodel.meas[wide].time,widths[wide],ps=2,co='aa'x
     oplot,!x.crange,[1,1]*0.1,co='aa00'x
     
     ; plot center
     delta = max(abs(centers-targetlinew[lineidx]))
     ; for MEGS-A sampling is 0.01 nm
     ; MEGS-B sampling is 0.02 nm
     cyr = targetlinew[lineidx] + ((3*dlambda)<delta>(2*dlambda))*[-1.,1.]
     plot,spmodel.meas.time, centers, xs=1,xr=xr,ys=1,yr=cyr,ps=-4, $
          xtit='Time (sec)',ytit='Center (nm)',$
          tit=strtrim(targetlinew[lineidx],2)+' nm'
     for ibin=-2,2 do oplot,!x.crange,ibin*dlambda + targetlinew[lineidx]*[1,1],linestyle=1
     oplot,!x.crange,[1,1]*targetlinew[lineidx],co='aa00'x,thick=2
     above=where(centers ge targetlinew[lineidx],n_above, comp=below, ncomp=n_below)
     if n_above gt 1 then oplot,spmodel.meas[above].time, centers[above], ps=6,co='aa'x
     if n_below gt 1 then oplot,spmodel.meas[below].time, centers[below], ps=2,co='ff0000'x
     
     ; plot offsets
     delta = max(abs(offsets))
     oyr = [-1.,1.]*delta
     plot,spmodel.meas.time, offsets, xs=1,xr=xr,yr=oyr,ys=1,ps=-4, $
          xtit='Time (sec)',ytit='Offsets (cps)',$
          tit=strtrim(targetlinew[lineidx],2)+' nm'
     low=where(offsets le 0.,n_low,comp=high,ncomp=n_high)
     if n_high gt 1 then oplot,spmodel.meas[high].time,offsets[high],ps=6,co='aa00'x
     if n_low gt 1 then oplot,spmodel.meas[low].time,offsets[low],ps=2,co='aa'x
     oplot,!x.crange,[0.,0.],co='aa00'x


     print,'median centers=',strtrim(median(centers),2),' widths=',strtrim(median(widths),2)
     ;stop                       ; use this stop to interrogate each line
        

     if !P.multi[0] eq 0 then stop
  endfor
  
  stop
  ; reset back to original plot settings
  !p.multi = pm_in
  !p.charsize= pc_in
  return
end

  
;+
; This procedure calculates a curve fit of each spectral wavelength vs
; time. The fitted function is a constant plus an exponential to match
; the spectral increase during the ascent. The physical meaning of the
; additive constant is the spectral irradiance at the top of the atmosphere
; where there is no atmospheric absorption. In practice, this absorption
; is expected to be less than 10% in most places near apogee (T+270), but
; could be around 20% in some places. 
;
; :Keywords:
;    spmodel: out, optional, type=structarr
;       This is the array of model results
;-
pro fit_megsa_36389_absorption, spmodel=spmodel

  !p.charsize=1.5 & !p.background='ffffff'x & !p.color=0
  @config36389 ; only needed for plots

  workingdir=file_dirname(routine_filepath()) ; in /code
  datadir = file_dirname(workingdir)+'/data/'

  radar = read_dat(datadir+'Radar_Data_36_389.dat')
  radartime_seconds = reform(radar[1,*])
  altitude_meters = reform(radar[9,*])
  
  restore,'rocket36389_megsb_irr.sav' ; spectra, spectra_cps, etc
  mb_spectra_cps = temporary(spectra_cps)
  mb_spectra = temporary(spectra)
  
  restore,'rocket36389_megsa_irr.sav' ; spectra, spectra_cps, etc
  
  ; try to fit an exponential?
  ; only use data prior to apogee of 299.34030 km (~T+283)
  ; roll started during index 247 (T+381) and finished just before 250 (T+411)
  ; index 256 (T+471) is likely partial with door closing

  ascending=where(spectra.time gt 120 and spectra.time lt 280)

  descending=where(spectra.time gt 410 and spectra.time lt 470)

  n_wave = n_elements(spectra[0].sp2)
  n_mbwave = n_elements(mb_spectra_cps[0].sp)
  sprec = { time:0., $
            sp1:fltarr(n_wave), fit1_result:fltarr(n_wave), $
            sp2:fltarr(n_wave), fit2_result:fltarr(n_wave), $
            mbsp:fltarr(n_mbwave), mbfit_result:fltarr(n_mbwave)  }
  spmodel = { predict_sp1:fltarr(n_wave), fit1_err:fltarr(n_wave), $
              predict_sp2:fltarr(n_wave), fit2_err:fltarr(n_wave), $
              wave1: spectra[0].w1, $
              wave2: spectra[0].w2, $
              wavemb: mb_spectra[0].w, $
              meas:replicate(sprec, n_elements(ascending)) }
  spmodel.meas.time = spectra[ascending].time
  spmodel.meas.sp1 = spectra_cps[ascending].sp1
  spmodel.meas.sp2 = spectra_cps[ascending].sp2
  spmodel.meas.mbsp = mb_spectra_cps[ascending].sp

  ; look for wavelength shift
  study_lambda, spmodel, /megsb
  stop
  study_lambda, spmodel, /slit1
  stop
  study_lambda, spmodel, /slit2
  
  ; assume the atmosphere density follows an equation similar to
  ; rho = rho0 * exp( -C * (altitude-a0) )
  ; where C has units of km^-1
  ; so rho/rho0 = exp( -C * (altitude-a0) )
  ; here we know the apogee is 293.98 from the config
  ; use 36.353 radar data to approximate alt at times
  ;sampletime_asc = spectra[ascending].time ; center time
  ; altitudes are km above sea level
  stop
  stop
  stop
  alt_asc = [180.140, 194.13, 207.131, 219.227, 230.387, 240.552, 249.99, 258.44, 265.902, 272.639, 278.367, 283.182, 287.04, 290.089, 292.246, 293.562]
  a0 = 1.26907 ; lowest value is last in descent at T+887 from 36.353, this is an upper limit

  ; functional form is f = A -B*exp(-Cx)

  ;widx = 3380                   ; a wavelength index into spectrum to fit
  for widx=200,3590 do begin
     ; slit2
     ;x = spectra_cps[ascending].time      ; seconds since launch
     x = alt_asc
     y = spectra_cps[ascending].sp2[widx] ; cps
     err = sqrt(abs(y)>1.)

     cpsfit2 = get_fit_params(x, y, err, perror, widx=widx, wave=spectra[0].w2[widx])

     spmodel.predict_sp2[widx] = cpsfit2[0]
     spmodel.meas.fit2_result[widx] = fitfunc(x, cpsfit2)

     spmodel.fit2_err[widx] = perror[0]
     
     ; slit 1
     ;x = spectra_cps[ascending].time   ; seconds since launch
     x = alt_asc
     y = spectra_cps[ascending].sp1[widx] ; cps
     err = sqrt(abs(y)>1.)

     cpsfit1 = get_fit_params(x, y, err, perror, widx=widx, wave=spectra[0].w2[widx])

     spmodel.predict_sp1[widx] = cpsfit1[0]
     spmodel.meas.fit1_result[widx] = fitfunc(x, cpsfit1)
     spmodel.fit1_err[widx] = perror[0]

     ;stop
  endfor

  refidx=[12,13,14,15]
  plot,spmodel.wave1,mean(spmodel.meas[refidx].sp1,dim=2),/ylog,yr=[1,1e4], $
       xtit='Wavelength (nm)', ytit='Slit 1 cps', $
       xrange=[5,20],xs=1, $
       tit='36.'+numberstr+' MEGS-A'
  oplot,spmodel.wave1,spmodel.predict_sp1,co='fe'x
  xyouts,8,10,'Measured mean'
  xyouts,8,20,'Extrapolation to zero absorption',color='fe'x
  stop
  plot,spmodel.wave1,mean(spmodel.meas[refidx].sp1,dim=2)/spmodel.predict_sp1,yr=[.7,1.1], $
       xtit='Wavelength (nm)', ytit='Extrap / Slit 1', $
       xrange=[5,20],xs=1, $
       tit='36.'+numberstr+' MEGS-A No atm prediction/meas mean'
  stop

  plot,spmodel.wave2,mean(spmodel.meas[refidx].sp2,dim=2),/ylog,yr=[1,1e5],$
       xtit='Wavelength (nm)', ytit='Slit 2 cps', $
       xrange=[10,40],xs=1, $
       tit='36.'+numberstr+' MEGS-A'
  oplot,spmodel.wave2,spmodel.predict_sp2,co='fe'x
  xyouts,18,10,'Measured mean'
  xyouts,18,20,'Extrapolation to zero absorption',color='fe'x
  stop
  plot,spmodel.wave2,mean(spmodel.meas[refidx].sp2,dim=2)/spmodel.predict_sp2,yr=[.7,1.1], $
       xtit='Wavelength (nm)', ytit='Extrap / Slit 2', $
       xrange=[10,40],xs=1, $
       tit='36.'+numberstr+' MEGS-A No atm prediction/meas mean'
  stop

  plot,spmodel.wave2,mean(spmodel.meas[refidx].sp2,dim=2)/spmodel.predict_sp2,yr=[.6,1.05], ps=-4,$
       xtit='Wavelength (nm)', ytit='Extrap / Slit 2', $
       xrange=[30,38],xs=1, $
       tit='36.'+numberstr+' MEGS-A No atm prediction/meas mean'
  stop

  ; plot 36-37 nm
  msp = mean(spmodel.meas[refidx].sp2,dim=2) ; mean measurement spectrum
  rsp = msp / spmodel.predict_sp2 ; mean measurement near apogee/model
  plot, spmodel.wave2, rsp, yr=[.6,1.1],/ys, ps=-4, xr=[36,37],/xs, $
        xtit='Wavelength (nm)', ytit='Ratio', $
        tit='Mean measurement/extrapolation'
  rp2 = spmodel.predict_sp2 / $
        max(spmodel.predict_sp2[where(spmodel.wave2 gt 36 and spmodel.wave2 lt 37)])
  rmsp = msp / max(msp[where(spmodel.wave2 gt 36 and spmodel.wave2 lt 37)])
  oplot,spmodel.wave2, rmsp*.1 + 1.,ps=-4,co='fe'x ; scale and offset
  xyouts,36.05,.7,'Mean scaled spectrum',co='fe'x
  xyouts,36.05, .67,'Measurement mean / extrapolated spectrum'
  stop
  
  plot,spmodel.wave2,spmodel.meas[14].sp2/spmodel.meas[15].sp2,xr=[36,37],ps=-4
  stop

  ;sampletime_des = spectra[descending].time ; center time
  ;alt_des = [216.788, 204.579, 191.404, 177.209, 162.225, 146.267]
  
  stop

  
  stop
  return
end
