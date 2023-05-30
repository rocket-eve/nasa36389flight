;+
; NAME:
;	megsb_sp_wave_cal
;
; PURPOSE:
;	Generate a wavelength scale using known lines in a solar spectrum
;
; CATEGORY:
;	EGS Level 1
;
; CALLING SEQUENCE:  
;	wavelength = megsb_sp_wave_cal( spectrum, err, status, file=file )
;
; INPUTS:
;	spectrum = solar spectrum (assumed fltarr(1024))
;	/file = optional specification of calibration file to read
;
; OUTPUTS:  
;	wavelength output is in nm
;	err = uncertainty [in nm] of wavelength scale to fitted line centers
;	status = status of result: 0 = OK, -1 = Invalid parameters,
;			-2 = Calibration file error, -6 = Bad calibration fit
;	/param = optional output for the 2nd order poly fit parameters
;
; COMMON BLOCKS:
;	megsb_sp_wave_com
;	------------
;		megsb_wave_set = TRUE if reference lines have been loaded
;		ref_wave = array of structure definition for reference lines
;			ref_wave.wave = parameter used here
;
; PROCEDURE:
;	1.  Check that parameters are valid
;	2.  Read calibration file if given in /file option of if not read yet
;	3.  For each line, find its center (assuming first guess wavelength scale)
;			A.  First search for uniform shift of spectrum
;			B.  Find location of maximum counts for each line
;			C.  Remove any common lines (lines blended too much)
;			D.  Calculate center of mass for line
;	4.  Fit wavelength scale to line center result
;			A.  First fit linear function
;			B.  Remove any wild points
;			C.  Finally fit 2nd order polynomial with cleaned data
;	5.  Return results
;
; MODIFICATION HISTORY:
;	7/17/1995	Hipook Brown	determine_egs_wave.pro  for METEOR SEE EGS
;	11/6/1999	Tom Woods		Modified for TIMED SEE EGS usage
;	02/04/2002	Don Woodraska	returns bad_cal_fit when shifting
;	causes k1 to exceed the length of temp_array
;	
;	FUTURE WORK:  fix slopes of 2nd order poly fit and only fit wavelength offset
;
; $Log: egs_wave_cal.pro,v $
; Revision 8.0  2005/06/15 18:50:13  see_sw
; commit of version 8.0
;
; Revision 8.0  2004/07/20 20:12:53  turkk
; commit of version 8.0
;
; Revision 7.0  2004/07/08 23:01:58  turkk
; commit of version 7.0
;
; Revision 6.0  2003/03/05 19:27:39  dlwoodra
; version 6 commit
;
; Revision 5.20  2002/09/06 23:18:01  see_sw
; commit of version 5.0
;
; Revision 4.0  2002/05/29 18:07:34  see_sw
; Release of version 4.0
;
; Revision 3.1  2002/02/06 18:58:59  see_sw
; added return if search range exceeds array bounds (bug fix)
;
; Revision 3.0  2002/02/01 18:54:35  see_sw
; version_3.0_commit
;
; Revision 1.6  2001/10/31 21:57:38  dlwoodra
; added return status of BAD_CAL_FIT when spectrum is too dark and cannot be fit.
;
; Revision 1.5  2001/10/30 19:21:17  dlwoodra
; added log file info
;
; Revision 1.4  2001/03/13 19:44:46  dlwoodra
; speed improvements only, vectorized xmin and xmax, removed IF statements from FOR loops
;
; Revision 1.3  2001/03/13 16:35:55  dlwoodra
; modified to avoid divide by zero for test data (search for machine_f.xmin)
;
; Revision 1.2  2001/02/13 00:53:01  dlwoodra
; After convolution added a check for channels above zero intensity. Returns reference wavelengths if less than 20 channels have positive intensity.
;
; Revision 1.1.1.1  2000/11/21 21:49:14  dlwoodra
; SEE Level 1 EGS Code Import
;
;
;idver='$Id: egs_wave_cal.pro,v 8.0 2005/06/15 18:50:13 see_sw Exp $'
;
;-

function megs_sp_wave_cal2, spectrum, err, status, file=file, param=param, megsb=megsb,debug=debug, slit1=slit1, order=order

common	megs_sp_wave_cal2, megsb_wave_set, ref_wave, cal_wave, last_file

if size(order,/type) eq 0 then order=3
; order is the order to pass into the polynomial fit, 2 means quadratic

;
;	Define CONSTANTS  (status values)
;
;@$see_code_hdr/see_defines.pro
;@$see_code_hdr/egs_defines.pro
@megsb_defines.pro

if keyword_set(debug) eq 1 then debug=1 else $
  DEBUG = 0                     ; set to non-zero to print DEBUG messages

;debug=1
;
;	setup MEGS-B wavelength scale initial guess (pre-calibration)
;
if keyword_set(megsb) then begin
;dwave = MEGSB_GRATING_D * MEGSB_ANODE_WIDTH / MEGSB_FOCAL_LEN
;    dwave = (97.702-40.193)/(1789-205) ; best guess using eyeball
    dwave = (97.702-49.9406)/(1813-500) ; best guess using eyeball 36.389
;
;pre_wave = 26.1 + dwave * findgen(2048)	      ; initial guess
    ; He guess
;    pre_wave = 31.679 + dwave * findgen(NUM_MEGSB_COLS) ;initial guess
;    param = [ 31.679, dwave ]
    ; Ne guess
    ;pre_wave = 35.2 + dwave * findgen(NUM_MEGSB_COLS) ;initial guess
    ;param = [ 35.2, dwave ]
    ;pre_wave = 32.677 + dwave * findgen(NUM_MEGSB_COLS) ;initial guess
    ;param = [ 32.677, dwave ]
    ;; Ne guess from centering in MOBI (2007066)
    ;pre_wave = 32.425 + dwave * findgen(NUM_MEGSB_COLS) ;initial guess
    ;param = [ 32.425, dwave ]
    ; need a better initial guess for 36.389
    tmpw = [49.9406, 58.4334, 97.702]
    tmppix = [501, 734, 1813]
    tmp = poly_fit(tmppix, tmpw,2)
    pre_wave = tmp[0] + tmp[1]*findgen(2048) + tmp[2]*(findgen(2048)^2)
    param = tmp

 endif else begin
    ; MEGS-A needs to be better defined

    ;rocket guesses (probably off-center)
    ;dwave = (30.3783-24.3027)/(1624-1288) ; best guess using eyeball
    ;pre_wave = 1.013 + dwave * findgen(NUM_MEGSB_COLS) ;initial guess
    ;param = [ 1.013, dwave ]
    ; rocket 36.258
    pre_wave = 4.22647 + (0.0134796 + findgen(NUM_MEGSB_COLS)*(1.58106e-06)) * findgen(NUM_MEGSB_COLS) ;initial guess
    dwave = (30.3783-17.107)/(1628.1273-867.11467) ; best guess using eyeball

    ; rocket 36.258
    ; 17.1 @ 868
    ; 18.06 @ 925
    ; 28.4 @ 1522
    ; 30.378 @ 1628
    ;hand_pix  = [868.,  925., 1522., 1628.] ; 36.358
    hand_pix  = [868.,  925., 1523., 1629.] ; 36.389
    hand_wave = [17.1, 18.06,  28.4, 30.378]
    tmp = poly_fit(hand_pix,hand_wave,2)
    print,'INFO: initial wavelength guess for slit 2 ',tmp
    pre_wave2 = tmp[0] + tmp[1]*findgen(2048) + tmp[2]*(findgen(2048)^2)
    ;stop

    ;mobi tank spectrum for flight EVE MEGS-A
    ;dwave = (30.3783-24.3027)/(1647-1311) ; best guess using eyeball
    ;pre_wave = 0.597 + dwave * findgen(NUM_MEGSB_COLS) ;initial guess
    ;param = [ 0.597, dwave ]

    ; FLIGHT DATA
    ;dwave = (30.3783-17.107)/(1648.55-891.28) ; best guess using eyeball
    ;pre_wave = 4.23879 + (0.0135944 + findgen(NUM_MEGSB_COLS)*(1.54202e-06)) * findgen(NUM_MEGSB_COLS) ;initial guess
;    pre_wave = 3.69372 + (0.0137110 + findgen(NUM_MEGSB_COLS)*(1.50176e-06)) * findgen(NUM_MEGSB_COLS) ;initial guess
    ; 2010/085
    ;pre_wave = 3.75975 + (0.0136065 + findgen(NUM_MEGSB_COLS)*(1.53837e-06)) * findgen(NUM_MEGSB_COLS) ;initial guess
    ; 2010/099
    ;pre_wave = 3.75719 + (0.01360003 + findgen(NUM_MEGSB_COLS)*(1.53856e-06)) * findgen(NUM_MEGSB_COLS) ;initial guess

    if keyword_set(slit1) then begin
       ;pre_wave = 4.40975 + (0.0133465 + findgen(NUM_MEGSB_COLS)*(1.66756e-06)) * findgen(NUM_MEGSB_COLS) ;initial guess
       ;good
       ;pre_wave = 4.39257 + (0.0134210 + findgen(NUM_MEGSB_COLS)*(1.60752e-06)) * findgen(NUM_MEGSB_COLS) ;initial guess
       ; test
       ;pre_wave = 4.40338 + (0.0133377 + findgen(NUM_MEGSB_COLS)*(1.70093e-06)) * findgen(NUM_MEGSB_COLS) ;initi
       ; FLIGHT GUESS
       ;pre_wave = 3.69372 + (0.0137110 + findgen(NUM_MEGSB_COLS)*(1.50176e-06)) * findgen(NUM_MEGSB_COLS) ;initial guess
       ; FLIGHT estimate from Harry's spectrum on 4/6/10
       ;pre_wave = 3.86799 + (0.0135148 + findgen(NUM_MEGSB_COLS)*(1.51956e-06)) * findgen(NUM_MEGSB_COLS) ;initial guess
       ; 36.258
       pre_wave = 4.44838 + (0.0134807 + findgen(NUM_MEGSB_COLS)*(1.53442e-06)) * findgen(NUM_MEGSB_COLS) ;initial guess

       ;; 36.389
       ;hand_pix  = [715., 857.] ; 36.389 TODO: revisit with known spectra
       ;hand_wave = [14.8377, 17.1]
       ;tmp = poly_fit(hand_pix,hand_wave,2)

    endif
    param = [ 0.597, dwave ]
endelse

err = 1.0				; in nm
status = BAD_PARAMS

;******************************************************************************
;
;	1.  Check that parameters are valid
;
if (n_params() lt 1) then begin
	print, 'Usage: wave = egs_wave_cal( spectrum, [ err, status, file=file, param=param ] )'
	return, pre_wave
endif

sp_size = size(spectrum)
if ( ( sp_size[0] ne 1 ) or ( sp_size[1] ne NUM_MEGSB_COLS ) ) then begin
	if (DEBUG ne 0) then print, 'ERROR: megsb_sp_wave_cal() requires "spectrum" to be '+strtrim(NUM_MEGSB_COLS,2)+' elements'
	return, pre_wave
endif

anode_nums = findgen(NUM_MEGSB_COLS)

;if ( (slit ne 0) or (slit ne 1) ) then begin
;	;  do nothing because not using "slit" information now anyway
;endif

;******************************************************************************
;
;	2.  Read calibration file if given in /file option or if not read yet
;		Return status=-2 if calibration file is not found.
;
read_the_file=0
if size(last_file,/type) eq 0 then last_file = ''
if size(file,/type) ne 0 then begin
    if file ne last_file then begin
        read_the_file=1
        last_file=file
    endif
endif
if (keyword_set(file) and read_the_file ne 0) then begin
	if (DEBUG ne 0) then print, 'Reading ' + file + ' ...'
	ref_wave = read_dat( file )
	loglun = get_log_lun()
    printf, loglun, ''
    printf, loglun, 'Number of REFERENCE LINE Records = ', $
      strtrim(n_elements(ref_wave),2)
	printf, loglun, 'MEGS_SP_WAVE_CAL: Read ' + file
endif else begin
	cr_size = size(megsb_wave_set)
	if (cr_size[cr_size[0]+1] eq 0) then megsb_wave_set = 0
	if (megsb_wave_set eq 0) then begin
        cal_dir = './' ; getenv( 'eve_cal_data' ) + '/' ; for SUN OS
		cal_file = 'megsb_ref_lines.dat'
		full_file = cal_dir + cal_file
		if (DEBUG ne 0) then print, 'Reading ' + full_file + ' ...'
		ref_wave = read_dat( full_file, /silent )
		loglun = get_log_lun()
        printf, loglun, ''
        printf, loglun, 'Number of REFERENCE LINE Records = ', $
          strtrim(n_elements(ref_wave),2)
		printf, loglun, 'MEGS_SP_WAVE_CAL: Read ' + full_file
	endif
endelse

if (n_elements(ref_wave) ge 2) then begin
	cal_wave = ref_wave[where(ref_wave.cal ne 0)]
endif else begin
	cal_wave = 0.
endelse

if (n_elements(cal_wave) lt 2) then begin
	status = BAD_FILE
	return, pre_wave
endif
megsb_wave_set = 1

status = OK_STATUS		; OK, passed all the initial checks

;******************************************************************************
;
;	3.  For each line, find its center (assuming initial wavelength scale)
;			A.  First search for uniform shift of spectrum
;			B.  Find location of maximum counts for each line
;			C.  Remove any common lines (lines blended too much)
;			D.  Calculate center of mass for line
;
   ; Define Wavelength Calibration Parameters
   max_anode_shift = 6    ; The Maximum Anodes the Spectrum Can Shift
   search_width    = 6    ; Half-Interval for Max Counts Search
   width           = 2.0  ; Half-Interval for Center of Mass Calculation
   max_yband       = 0.4  ; Defines Error Tolerance for Centers of Mass in nm

   ; Predict Reference Anode for Reference Wavelengths
   n_ref = n_elements(cal_wave)
   ref_anodes = fltarr( n_ref ) ; same size as reference lines
   ref_anodes = round( interpol(findgen(NUM_MEGSB_COLS), pre_wave, cal_wave.wave) )

   ;Create a Reference Spectrum
   ref_spectrum = fltarr(NUM_MEGSB_COLS) ; some big array
   ref_spectrum[ref_anodes] = 1.0
   ref_offset = min(ref_anodes)
   ref_spectrum = ref_spectrum[ min(ref_anodes) : max(ref_anodes) ]
   
   ;Convolve Reference Spectrum with Counts
   temp_array = convol(abs(spectrum),ref_spectrum)
   if (n_elements(where(temp_array ne 0)) lt 20) then begin
       print,'EGS_WAVE_CAL: cannot fit wavelengths to spectrum, too dark'
       status=BAD_CAL_FIT
       return,pre_wave
   endif
   temp_array = temp_array[where(temp_array ne 0)]

   ;Create Filter Window
   window = findgen(n_elements(temp_array))
   center = ref_offset  ; centered around theoretical
   window = -(window - center)^2 ; window is an inverse quadratic
   window = (window - min(window))
   temp_array = temp_array*window

   ;Find the Max of the Convolution & Create New Ref_Anodes
   theMax = max( temp_array, anode_shift )
   anode_shift = anode_shift - ref_offset
   if (DEBUG ne 0) then begin
   		print, 'egs_wave_cal() found uniform shift of ', $
   			anode_shift,  ' anodes.'
   	endif
   	
   ;If Anode Shift is Too Large, then restrict search for peak
   if abs(anode_shift) gt 2*search_width then begin
   		k1 = (ref_offset - 2*search_width) > 0 ;no less than zero
   		;if k1 lt 0 then k1 = 0
   		k2 = (ref_offset + 2*search_width) < (n_elements(temp_array) - 1)
   		;if k2 ge n_elements(temp_array) then k2 = n_elements(temp_array)-1

        ;
        ;DLW 2-4-02
        ;
        if k1 gt n_elements(temp_array)-1 then begin
            print,'EGS_WAVE_CAL: error, wavelength peak is out of range'
            status=BAD_CAL_FIT
            return,pre_wave
        endif
   		theMax = max( temp_array[k1:k2], anode_shift )
   		anode_shift = anode_shift + k1 - ref_offset
		if (DEBUG ne 0) then begin
	   		print, 'egs_wave_cal() using a shift of ', anode_shift, ' anodes.'
	   		;stop, 'STOPPED to check out anode_shift and temp_array problems...'
	   	endif
   endif

	if (DEBUG ne 0) then begin
		plot, temp_array, title = 'Anode Shift = '+strtrim(anode_shift,2)
		oplot, [anode_shift+ref_offset,anode_shift+ref_offset], $
			[!cymin, !cymax], line=2
   		ans = ''
   		read, 'Continue ? ', ans
   	endif
   
   ref_anodes_temp = ref_anodes + anode_shift

   ;Find Anodes Where Nearest Max Count Occurs
   max_anodes = ref_anodes_temp ; Array of Anodes Where Max Counts Occur

   ;Define Search Area
   x_min = (ref_anodes_temp - search_width) > 0 ;no less than zero
   x_max = (ref_anodes_temp + search_width) < NUM_MEGSB_COLS
   for i = 0, n_ref-1 do begin
      ;anode = ref_anodes_temp[i]

      ;Define Search Area
      ;if (x_min lt 0) then x_min = 0
      ;if (x_max ge NUM_MEGSB_COLS) then x_max = NUM_MEGSB_COLS
      
      ;Get Nearest Maximum Count
      theMax = max(spectrum(x_min[i]:x_max[i]), pick )
      max_anodes[i] = x_min[i] + pick
   endfor

   ;See If Any Max_Anodes Were Picked Twice
   ;double_check_dn = fltarr(n_ref) + 10 ; initialize to 10
   ;double_check_up = double_check_dn
   ;for i = 0, n_ref-2 do begin
   ;   double_check_dn[i+1] = abs(max_anodes[i+1] - max_anodes[i]) 
   ;   double_check_up[i]   = abs(max_anodes[i] - max_anodes[i+1])
   ;endfor

   ; the following provides no noticeable speed improvement
   double_check_dn = abs(shift(max_anodes,-1) - max_anodes)
   double_check_up = abs(max_anodes - shift(max_anodes,-1))
   double_check_dn[0] = 10       ;revert element 0 to initialization
   double_check_up[n_ref-1] = 10 ;revert last element to initialization

   err_width = search_width
   d_pick = where( (double_check_dn le err_width) or (double_check_up le err_width) )

   ;If Max_Anode Picked Twice, Then Recalculate Max Anodes
   if d_pick[0] ne -1 then begin
      n_dpick = n_elements(d_pick)
	  if (DEBUG ne 0) then begin
	  	print, 'egs_wave_cal: Double peaks found for ', n_dpick, ' lines: '
	  	print, cal_wave[d_pick].wave
	  endif
	  
      ; Reduce Search Width to Insure no Anodes Will Be Picked Twice
      search_width_2 = max(abs(max_anodes[d_pick]-ref_anodes_temp[d_pick])) / 2.

	  ;Define Search Area
	  x_min = (ref_anodes_temp - search_width_2) > 0 ;no less than zero
	  x_max = (ref_anodes_temp + search_width_2) < NUM_MEGSB_COLS ;at most 1024
      for j = 0, n_dpick-1 do begin

		i = d_pick[j]
		;anode = ref_anodes_temp[i]
		
		;Define Search Area
		;x_min = anode - search_width_2
		;x_max = anode + search_width_2
		;if (x_min lt 0) then x_min = 0
		;if (x_max ge NUM_MEGSB_COLS) then x_max = NUM_MEGSB_COLS
		
		;Get Nearest Maximum Count
		theMax = max(spectrum(x_min[i]:x_max[i]), pick )
		max_anodes[i] = x_min[i] + pick
      endfor

   endif

   ;Calculate Centers of Mass for Each Maximum Anode

   CM_anodes = float(max_anodes) 		; array of centers of mass
   CM_values = spectrum[max_anodes]		; array of masses
   CM_anodes_2 = CM_anodes

   ;First This Calculates the Center of Mass Around the Max_Anodes, Then
   ;it Calculates the Center of Mass Around the Center of Mass 2 more times.
   for j = 1, 3 do begin

      ;Search Area
      low  = round(CM_anodes - width) > 0 ;no less than zero
      high = round(CM_anodes + width) < NUM_MEGSB_COLS ;no more than 1024
      for i = 0, n_ref-1 do begin
         ;Search Area
         ;low  = CM_anodes[i] - width
         ;high = CM_anodes[i] + width
         ;Do Not Calculate Beyond Limits
         ;if (low ge 0) and (high le NUM_MEGSB_COLS) then begin
         ;   low  = round(low)
         ;   high = round(high)
         ;   ;Center of Mass Calculation
            tot_spec=total(spectrum[low[i]:high[i]])
            ; compare tot_spec to smallest positive floating point number
            ; to avoid a divide by zero (NaN)
            if (tot_spec ge machine_f.xmin) then CM_anodes_2[i] = $
                total(spectrum[low[i]:high[i]]*anode_nums[low[i]:high[i]]) $
                  / tot_spec
            ;stop,'egs_wave_cal: stopped'
         ;endif
      endfor

      ;If the Center of Mass Around the Center of Mass Differs Less Than 
      ;Err_Width Anode then Keep the New Center of Mass
      diff = abs(CM_anodes_2 - CM_anodes)
      err_width = 1
      good_pick = where(diff le err_width)
      bad_pick  = where(diff gt err_width)

      ;Return Error If no Good_Picks
      if n_elements(good_pick) lt 2 then begin
         if (DEBUG ne 0) then print, 'ERROR: egs_wave_cal() failed for Center of Mass calculation !'
         status = BAD_CAL_FIT
         return, pre_wave
      endif

      ;Keep Good CM_anodes
      if good_pick(0) ne -1 then CM_anodes[good_pick] = CM_anodes_2[good_pick]
      ;Mark Bad CM_anodes by -1
      if bad_pick(0)  ne -1 then CM_anodes[bad_pick] = -1
	  
   endfor

if (DEBUG ne 0) and (bad_pick[0] ne -1) then $
	print, 'egs_wave_cal: Number of bad Center of Mass values = ', n_elements(bad_pick)
	
;******************************************************************************
;
;	4.  Fit wavelength scale to line center result
;			A.  First fit linear function
;			B.  Remove any wild points
;			C.  Finally fit 2nd order polynomial with cleaned data
;
   ;Get New Calibration Coefficients

   ;Do Initial Linear Fit
   pick = where(CM_anodes ge 0) ; Remove Bad Data Points

   c1 = poly_fit(CM_anodes[pick], cal_wave[pick].wave, 1, Yfit1, $
                   Yband1, fit1_err )
   ; use the pre_wave instead 7-24-13 DLW
   yfit1 = pre_wave[round(cm_anodes[pick])]
;   stop

   ;Discard Bad Data Points, Where Yband Less Than Max_Yband
   CM_anodes1 = CM_anodes[pick]
   pick2 = where(abs(Yfit1 - cal_wave[pick].wave) ge max_yband)
   if pick2[0] ne -1 then CM_anodes[pick[pick2]] = -2
   pick1 = pick
   pick = where(CM_anodes ge 0)

   ;Report Bad CM_anodes
   if (DEBUG ne 0) and (pick2[0] ne -1) then $
      print, 'egs_wave_cal: Number of wild points = ', n_elements(pick2)

	;param = c1
        ;new_wave = c1[0] + c1[1] * anode_nums
   new_wave = pre_Wave 
     
   if (DEBUG ne 0) then begin
      loadct, 12
      plot, CM_anodes1, cal_wave[pick1].wave, xtitle='Anode', $
            ytitle='Wavelength (nm)', psym=4,ys=1
      oplot, CM_anodes1, Yfit1, line=2, psym=0, color=150
      if (n_elements(pick2) gt 1) then $
         oplot, CM_anodes1[pick2], cal_wave[pick1[pick2]].wave, $
                psym=5, color=150 ; pick2 are wild points
   endif
	
	;Return Error If Not Enough Picks
   if n_elements(pick) le 3 then begin
      if (DEBUG ne 0) then print, 'WARNING: egs_wave_cal() only did a linear fit !'
      status = BAD_CAL_FIT
      stop
      err = fit1_err
      return, new_wave
   endif

   ;Get the Calibration Coefficients From sine wave fit
   new_coeff = poly_fit(CM_anodes[pick], cal_wave[pick].wave, order, Yfit=yfit, $
                   Yband=yband, sigma=fit_err )
   start=[0.d, 1100.,.00001]
   rerr=fltarr(n_elements(pick)) + 1.
   dlw_new_coeff = mpfitfun('sine_fun', cm_anodes[pick], cal_wave[pick].wave,rerr, start, niter=2000 )

   start=[0.d, 1100.,.00001,.00001]
   dlw2_new_coeff = mpfitfun('sine2_fun', cm_anodes[pick], cal_wave[pick].wave,rerr, start,niter=2000 )
   if DEBUG ne 0 then begin
      oplot,CM_anodes1,sine_fun(cm_anodes1,dlw_new_coeff),ps=-4,co='fe'x
      oplot,CM_anodes1,sine2_fun(cm_anodes1,dlw2_new_coeff),ps=-2,co='fe00'x
      stop
   endif
   ;Get the Calibration Coefficients From Polynomial Fit
   ;new_coeff = poly_fit(CM_anodes[pick], cal_wave[pick].wave, 2, Yfit, $
   ;                Yband, fit_err )

   ;Get the Wavelengths for Each Anode Using New Calibration Coeff
   wavelengths = new_coeff[0]
   for i = 1, n_elements(new_coeff) - 1 do begin
      wavelengths = wavelengths + (float(anode_nums)^i)*new_coeff[i]
   endfor

	; save error (fit 1-sigma uncertainty) in nm
	param = new_coeff
	err = fit_err
	
	if (DEBUG ne 0) then begin
		print, ' '
		print, 'Linear fit parameters = ', c1[0], c1[1]
		print, '            fit error = ', fit1_err, ' nm'
		print, ' '
		;print, '2nd Order fit params  = ', new_coeff[0],  new_coeff[1],  new_coeff[2]
		;print, '            fit error = ', fit_err, ' nm'
		print, '3rd Order fit params  = ', new_coeff[0],  new_coeff[1],  new_coeff[2], new_coeff[3]
		print, '            fit error = ', fit_err, ' nm'
		print, ' '
		oplot, CM_anodes[pick], Yfit
		if (n_elements(pick2) gt 1) then $
			oplot, CM_anodes1[pick2], cal_wave[pick1[pick2]].wave, $
				psym=5, color=150
		ans = ''
		read, 'Hit ENTER key when ready for next plot ', ans
	endif
	
;******************************************************************************
;
;		Plot up results if in DEBUG mode
;
if (DEBUG ne 0) then begin
   ;Plot_EGS
   loadct, 12
   xrange = [min(wavelengths),max(wavelengths)]
   loop = 0
PLOT_R:
   plot, wavelengths, spectrum, xrange = xrange, psym=10

   ;Plot Reference Lines $ Max Lines $ Centers Mass
   x1 = wavelengths[round(ref_anodes_temp)]
   y1 = spectrum[ref_anodes_temp] * 0.9
   x2 = wavelengths[round(max_anodes)]
   y2 = spectrum[max_anodes]
   wgd = where( CM_anodes ge 0 )
   x3 = interpol( wavelengths, anode_nums, CM_anodes[wgd] )
   y3 = y2[wgd]/2
   oplot, x1, y1, color = 85,  psym = 5
   oplot, x2, y2, color = 120, psym = 6
   oplot, x3,  y3, color = 185, psym = 4

   oplot, sine_fun(findgen(2048),dlw_new_coeff), spectrum,co='fe'x
   oplot, sine2_fun(findgen(2048),dlw2_new_coeff), spectrum,co='aa00'x

   if loop then begin
      xyouts, x1, y1, alignment = -0.5, $
              strtrim(x1,2), color = 85
      xyouts, x2, y2, alignment = -0.5, $
              strtrim(x2,2), color = 120
      xyouts, x3, y3,  alignment = -0.5, $
              strtrim(x3, 2), color = 185
      y1 = !cymin
      y2 = !cymax
      for i = 0, n_ref-1 do begin
         x1 = wavelengths[round(ref_anodes_temp[i])]
         x2 = wavelengths[round(max_anodes[i])]
;         x3 = wavelengths[round(CM_anodes[i])]
;         oplot, [x1, x1], [y1, y2], color = 85 , linestyle = 2
;         oplot, [x2, x2], [y1, y2], color = 120, linestyle = 4
         oplot, [x2, x2], [y1, y2], color = 185, linestyle = 1
      endfor
   endif

   ;See Zoomed Views
   r_number = -1
   read, r_number, prompt = 'Enter Reference Wavelength Index (-1 to exit) : '
   if (r_number ge 0) and (r_number lt n_ref) then begin
      x = wavelengths[round(max_anodes[r_number])]
      xrange = [x - 3*search_width*dwave, x + 3*search_width*dwave]
      loop = 1
      goto, PLOT_R
   endif
endif

;******************************************************************************
;
;	5.  Return results
;
; if (DEBUG ne 0) then stop, 'STOPPED at end of egs_wave_cal() ...'
return, wavelengths

end
