function get_36389_sensitivity, megsa1=megsa1, megsa2=megsa2, megsb=megsb

  datapath = '../data/'

  if keyword_set(megsa1) then begin
     ; sensitivity_MEGSA1_second_order_2017.sav has two additional vectors
     ; spec_corrected1: This vector contains spec_uncorrected corrected for second order using 2009 data
     ; spec_corrected2: This vector contains spec_uncorrected corrected for second order using 2010 data

     ;restore,'data/sensitivity_MEGSA1_second_order_2013.sav' ; sensitivity change from Brian 7/17/19 updated 7/19/19
     ; the above saveset has too low of values at 17.1 nm for both slits
     
     print,'WARNING: USING 36336(2018) MA2 SENSITIVITY'
     restore,datapath+'36336sensitivity_MEGSA1_second_order_2017.sav' ; sensitivity from Brian 7/17/19

     ; Brian reverses some sensitivity maps and not others
     ; zero is used to fill so use one row to find where the sensitivity
     ; data is located to figure out if reversing is needed
     if total(sens_corrected1[0:1023,700]) gt total(sens_corrected1[1024:2047,700]) then begin
        sens_corrected1 = reverse(sens_corrected1)
        sens_uncorrected = reverse(sens_uncorrected)
        sens_corrected2 = reverse(sens_corrected2)
     endif

     newsens = sens_corrected1  ; select this one, but replace bad parts
     
     ; remove long wavelength bad stuff in slit 1 beyond column 1130
     newsens[1130:*,*] = sens_uncorrected[1130:*,*] ; higher order correction fails in longer wavelengths
     ; that part is contaminated with higher orders at SURF
     
     ; fill missing data with a line fit in cross-dispersion direction
     for i=0,2047 do begin
        plot,newsens[i,*],/ylog,yr=[1e-8,.01],tit='MA1 col#'+strtrim(i,2),ytit='Sensitivity',xtit='Row'
        gd=where(newsens[i,*] gt 1.1e-8 and $
                 newsens[i,*] lt 0.01, n_gd, comp=comp)
        if n_elements(comp) gt 2 and n_gd gt 2 then begin
           ; both good and bad arrays exist
           coef=ladfit(gd,newsens[i,gd])
           ;if i eq 1130 then stop
           line = coef[0] + coef[1]*dindgen(2048)
           newsens[i,comp] = line[comp] ; fill vertically
           rat = newsens[i,*] / line ; ratio to linefit
           bad=where((rat gt 10) or (rat lt .1),n_bad,comp=good) ; rat should be close to 1.0
           if n_bad gt 0 then newsens[i,bad] = interpol(newsens[i,good],good,bad)
           oplot,newsens[i,*],co='fe'x ; corrected/filled
           oplot,sens_uncorrected[i,*],co='fe0000'x ; no higher order correction
           ;if i eq 1129 then stop
           
        endif
     endfor
     
     new1 = (newsens < .01) > (1e-8) ; 2017 380Mev cal corrected for 2nd order using 2009 data
     ;new1 is in wavelength order now 8/19/20
     ;ma1 side is new1[*,512:1023]
     
     return, new1
  endif

  if keyword_set(megsa2) then begin
     
     ; the file data/sensitivity_MEGSA2_second_order_v2_2013.sav is bad

     ; the file data/sensitivity_MEGSA2_fw0_380MeV_2013.sav has these things
     ;ALPHA           FLOAT     =       0.00000
     ;BEAMENERGY      INT       =      380
     ;BETA            FLOAT     =       0.00000
     ;FILTER          INT       =        0
     ;INSTRUMENT      STRING    = 'MEGSA2'
     ;SENSITIVITY     FLOAT     = Array[2048, 1024]
     ;SENSITIVITY_ERROR
     ;           FLOAT     = Array[2048, 1024]
     ; The sensitivity_error is all zeroes.
     ; Wavelengths are backwards (normal), long wave is at 0
     ; 
     
     
     ; REVERTING TO 36336 sensitivity map for A2
     ; sensitivity_MEGSA2_second_order_2017.sav has three additional vectors
     ; spec_corrected1: This vector contains spec_uncorrected corrected for second order using 2009 data
     ; spec_corrected2: This vector contains spec_uncorrected corrected for second order using 183MeV data from 2009
     ; spec_corrected3: This vector contains spec_uncorrected corrected for second order using 285MeV data from 2009

     print,'WARNING: USING 36336(2018) MA2 SENSITIVITY'
     restore,datapath+'36336sensitivity_MEGSA2_second_order_2017.sav' ; sensitivity from Brian 7/17/19
     ; MEGS-A2 save file from 36336 contains these items
     ;COMMENTS        STRING    = Array[4]
     ;SENS_CORRECTED1 FLOAT     = Array[2048, 1024] ; 
     ;SENS_CORRECTED2 FLOAT     = Array[2048, 1024]
     ;SENS_CORRECTED3 FLOAT     = Array[2048, 1024]
     ;SENS_UNCORRECTED
     ;                FLOAT     = Array[2048, 1024]
     ;WAVE            FLOAT     = Array[3600]

     ;restore,'data/sensitivity_MEGSA2_second_order_2013_2.sav' ; sensitivity from Brian 7/17/19

     ; Brian reverses some sensitivity maps and not others
     ; zero is used to fill so use one row to find where the sensitivity
     ; data is located to figure out if reversing is needed
     if total(sens_corrected1[500:600,200]) lt total(sens_corrected1[1500:1600,200]) then begin
        ; reverse 8/19/20
        sens_corrected1 = reverse(sens_corrected1)
        sens_uncorrected = reverse(sens_uncorrected)
        sens_corrected2 = reverse(sens_corrected2)
        sens_corrected3 = reverse(sens_corrected3)
     endif
     
     newsens = sens_corrected1 ; select this one, and may need to replace parts

     ; fill missing data with a line fit
     for i=0,2047 do begin
        plot,newsens[i,*],/ylog,yr=[1e-8,.01],tit='MA2 col#'+strtrim(i,2),ytit='Sensitivity',xtit='Row'
        gd=where(newsens[i,*] gt 1.1e-8 and $
                 newsens[i,*] lt 0.01, n_gd, comp=comp)
        if n_elements(comp) gt 2 and n_gd gt 2 then begin
           ; both good and bad arrays exist
           coef=ladfit(gd,newsens[i,gd])
           line = coef[0] + coef[1]*dindgen(2048)
           newsens[i,comp] = line[comp] ; fill vertically
           rat = newsens[i,*] / line ; ratio to linefit
           bad=where((rat gt 10) or (rat lt .1),n_bad,comp=good) ; rat should be close to 1.0
           if n_bad gt 0 then newsens[i,bad] = interpol(newsens[i,good],good,bad) ; fill vertically
           oplot,newsens[i,*],co='fe'x ; corrected/filled
           oplot,sens_uncorrected[i,*],co='fe0000'x ; no higher order correction

        endif
     endfor

     ; replace slit 1 area with a repeated reference value from the slit 2 side
     ref = mean(newsens[*,300:400],dim=2) ; mean of the good part of slit 2
     for i=512,1023 do newsens[*,i] = ref ; (not measured at SURF)
     
     new1 = (newsens < .01) > (1e-8)
     
     return,new1
  endif

  if keyword_set(megsb) then begin
     ; from Brian 7/19/19 v38 ***
     sensfile=datapath+'36336sensitivity_MEGSB_fw0_380MeV_2017.sav' ; scaled to A2 in overlap
     ; alpha, beta, filter, instrument, sensitivity[2048,1024], sensitivity_error[2048,1024]

     restore,sensfile
     return,sensitivity
     
  endif

  return,-1
end
