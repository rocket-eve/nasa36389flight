;
;   extract_megs_lines.pro
;
;   Make MEGS spectrum from MEGS image (using megs_makespectrum.pro)
;	and then process to find lines - plot results
;
;   8/14/06
;   Tom Woods
;
;   INPUT:
;	  w	         wavelength array
;	  sp		 spectrum array
;	  channel	 Options are  'B', 'A', 'A1', 'A2', 'SAM'
;	  width		 Option to specify width for line search (default is 400)
;     filename   Option to read MEGS FITS file (using readmegs.pro)
;     /debug     Option to display messages
;	  /save		 Option to save line list
;
;   OUTPUT:
;     Function returns spectrum
;
;	  linelist	List of found lines (line center and line peak)
;
;function  extract_megs_lines, w, sp, channel, width=width, filename=filename, debug=debug, $
;						linelist=linelist, save=save

function  extract_megs_lines, sp, width=width, filename=filename, debug=debug, $
						linelist=linelist, save=save

if (n_params() lt 1) and (not keyword_set(filename)) then begin
  print, 'USAGE: spectrum = megs_lines( image, channel [, filename=filename, /debug ] )'
  filename = ''
  print, ' '
endif 

;if keyword_set(filename) or ((not keyword_set(filename)) and (size(filename,/type) eq 7)) then begin
;  if keyword_set(debug) then image = readmegs( filename, /debug )  $
;  else image = readmegs( filename, /debug )    ; ignores "image" input
;  if (strlen(filename) gt 5) then begin
;    spos = strpos(filename,'/',/reverse_search)
;    if (spos ge 0) then sfile = strmid(filename,spos+1,strlen(filename)-spos)
;  endif
;endif else sfile = ''

;if n_elements(image) lt (2048L * 1024L) then begin
;  print, 'MEGS_LINES:  error in "image" size'
;  return, -1L
;endif

;if n_params() lt 2 then begin
;  channel = ' '
;  read, 'Enter channel type [B, A, A1, A2, SAM] ? ', channel
;endif

;ch = strupcase(strmid(channel,0,2))
;ch1 = strmid(ch,0,1)

;if (strlen(sfile) lt 1) then sfile = 'MEGS-'+ch

;;
;;	make spectrum now (using megs_makespectrum.pro)
;;
;if keyword_set(debug) then begin
;  spectrum = megs_makespectrum( image, ch, /background, /debug )
;endif else begin
;  spectrum = megs_makespectrum( image, ch, /background, /debug )
;endelse
spectrum=sp

dbins=n_elements(sp) ; 2048 with raw vertical img sum

if not keyword_set(width) then width = dbins
nplots = dbins / width
xbins = indgen(dbins)
linecenter = fltarr(width)
linepeak = fltarr(width)
linewidth = fltarr(width)
line1 = fltarr(width)
line2 = fltarr(width)
ans = ' '

linelist = fltarr(3,dbins)
linecnt = 0L

;LINE_WIDTH = round(5L * float(dbins)/2048L)
LINE_WIDTH = round(5L * float(dbins)/2048L)
LINE_MIN_HEIGHT = 300L


;  plot setup
;setplot
csize = 1.5

for k=0L,nplots-1 do begin
stop,'something is wrong in this code that makes it run away'
  x1 = k*width
  x2 = x1 + width
  cnt = 0L
  wrange = where( (xbins ge (x1-LINE_WIDTH*2)) and (xbins lt (x2+LINE_WIDTH*2)) )
  tsp = float(spectrum[wrange])
  tx = float(xbins[wrange])
  nsp = n_elements(tsp)
  tf = intarr(nsp) + 1
  tf[0] = 0
  tf[nsp-1] = 0
  ;  define MIN as lower quartel
  tmin = median(tsp)
  wgd1 = where( tsp lt tmin )
  tmin = median(tsp[wgd1])
  wgd1 = where( tsp lt tmin )  
  tmin = median(tsp[wgd1])
  if (tmin lt 0) then tmin=0
  tminline = tmin + LINE_MIN_HEIGHT
  wbad = where( tsp lt (tminline), nbad )
  if (nbad gt 0) then tf[wbad] = 0
  wgood = where ( (tsp ge (tminline)) and (tf ge 1), ngood )
  if (ngood gt 1) then begin
    ;
    ; find the lines and save the line peak and line center (center of mass)
    ;
    while (ngood gt 0) do begin
      tmax = max( tsp[wgood], wmax )
      imax = wgood[wmax]
      i1 = wgood[wmax]
      i2 = wgood[wmax]
      i1r = i1 - LINE_WIDTH*2
      if (i1r lt 1) then i1r = 1
	  i2r = i2 + LINE_WIDTH*2
	  if (i2r ge (nsp-1)) then i2r = nsp-2
	  for i=i1,i1r,-1 do begin
	    if (tsp[i] lt tminline) or (tf[i] lt 1) then begin
	      break		; break from the "for" loop
	    endif
	  endfor
	  i1 = i
	  for i=i2,i2r,1 do begin
	    if (tsp[i] lt tminline) or (tf[i] lt 1) then begin
	      break		; break from the "for" loop
	    endif
	  endfor
	  i2 = i
	  tf[i1:i2] = 0		; mark as being used already
	  if (i1 eq imax) then i1 = imax - 1
	  if (i2 eq imax) then i2 = imax + 1
	  linepeak[cnt] = tmax
	  linecenter[cnt] = total( tx[i1:i2] * tsp[i1:i2] ) / total( tsp[i1:i2] )
	  line1[cnt] = i1
	  line2[cnt] = i2

      ;hack to find width
      hw_thresh=1L
      lo=(imax-hw_thresh)>0L
      hi=(imax+hw_thresh)<(n_elements(tx)-2L)
      while tsp[lo] gt tmax/2. do lo=(lo-1L)>0 ;left edge
      while tsp[hi] gt tmax/2. do hi=(hi+1)<(n_elements(tx)-1) ;right edge
      leftedge=interpol(tx[lo:imax],tsp[lo:imax],tmax/2.d0)
      rightedge=interpol(tx[imax:hi],tsp[imax:hi],tmax/2.d0)
      linewidth[cnt]=rightedge*1.d0-leftedge
;      if linewidth[cnt] gt 10 or linewidth[cnt] lt 0 then linewidth[cnt]=0

	  cnt = cnt + 1L
      wgood = where ( (tsp ge (tmin+LINE_MIN_HEIGHT)) and (tf ge 1), ngood )
    endwhile
    ;
    ;	save the line results
    ;
    if (cnt gt 0) then begin
      linelist[0,linecnt:linecnt+cnt-1] = linecenter[0:cnt-1]
      linelist[1,linecnt:linecnt+cnt-1] = linepeak[0:cnt-1]
      linelist[2,linecnt:linecnt+cnt-1] = linewidth[0:cnt-1]
      linecnt = linecnt + cnt
    endif
    
    ;
    ; plot the results
    ;
    if keyword_set(debug) gt 0 then begin
        yr = [ 0., (max(tsp) + LINE_MIN_HEIGHT) * 1.2 ]
        plot, tx, tsp, xrange=[x1,x2], xs=1, psym=10, yrange=yr, ys=1, $
              xtitle='Pixel', ytitle='Signal (DN)' ;, title=sfile
        for i=0,cnt-1 do begin
            if (linecenter[i] gt x1) and (linecenter[i] lt x2) then begin
                oplot, linecenter[i]*[1.,1], linepeak[i]*1.02+[0.,LINE_MIN_HEIGHT], thick=1
                xyouts, linecenter[i], linepeak[i]*1.04 + LINE_MIN_HEIGHT, $
                        strtrim(string(linecenter[i],format='(F7.1)'),2), $
                        orient=90.,charsize=csize
            endif
        endfor
        if (k lt (nplots-1)) then read, 'Next plot ? ', ans
    endif
endif
endfor

if (linecnt lt 1) then linelist = -1L else begin
  linelist = linelist[*,0:linecnt-1]

  ;filter out lines that have bad widths
  tmp=where(linelist[2,*] gt 0,n_tmp)
  linelist = linelist[*,tmp]

  linesort = sort(reform(linelist[0,*]))
  temp = linelist
  linelist[0,*] = temp[0,linesort]
  linelist[1,*] = temp[1,linesort]
  linelist[2,*] = temp[2,linesort]
;  if keyword_set(save) and keyword_set(filename) then begin
;    savefile = filename + '.linelist'
;    if (keyword_set(debug)) then print, 'Saving line list in ', savefile
;    write_dat, linelist, file=savefile, $
;    		lintext='Line list from megs_lines.pro', $
;    		coltext='Line Center (pixel),  Line Peak (DN)'
;  endif
endelse
;stop
return, linelist ;spectrum
end
