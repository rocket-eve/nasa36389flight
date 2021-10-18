pro plot_mb_lines, wave, sp, lines

for i=0L,n_elements(lines)-1 do begin
   ;if 1 eq 1 or lines[i].cal ne 0 then begin
   if lines[i].cal ne 0 then begin
      xrange = lines[i].wave+[-1,1]
      tmp=min(abs(wave-lines[i].wave),ii)
      yrange=[100,max(sp[ii])*10]
      plot,wave,sp,xr=xrange,xs=1,/ylog,yr=yrange,ps=10
      idx=where(lines.wave gt xrange[0] and lines.wave lt xrange[1],n_idx)
      for j=0L,n_idx-1 do begin
         if lines[idx[j]].cal ne 0 then begin
            co='f0'x 
            li=2
            factor=1.5
         endif else begin
            co='ff'x
            li=1
            factor=.1
         endelse
         tmp=min(abs(wave-lines[idx[j]].wave),ii)
         oplot,[1,1]*lines[idx[j]].wave, [1,sp[ii]],lines=li,co=co
         xyouts,lines[idx[j]].wave,sp[ii]*factor,strtrim(string(lines[idx[j]].wave,form='(f8.3)'),2)+' '+lines[idx[j]].name,orient=90
      endfor
      stop
   endif
endfor

;for k=32,106,8 do begin
;   xrange=[k-1,k+8]
for k=32,106,4 do begin
   xrange=[k-1,k+4]

for i=0L,n_elements(lines)-1 do begin
   if lines[i].cal ne 0 then begin
;      xrange = lines[i].wave+[-1,1]*3
      tmp=min(abs(wave-lines[i].wave),ii)
;      yrange=[100,max(sp[ii])*10]
      yrange=[100,1e7]
      plot,wave,sp,xr=xrange,xs=1,/ylog,yr=yrange,ps=10
      idx=where(lines.wave gt xrange[0] and lines.wave lt xrange[1],n_idx)
      for j=0L,n_idx-1 do begin
         ion=(strsplit(lines[idx[j]].name,' ',/extract))[0]
         ;stop
         if lines[idx[j]].cal ne 0 or ion eq 'H' then begin
            co='f0'x 
            li=2
            factor=1.5
            name=strtrim(string(lines[idx[j]].wave,form='(f8.3)'),2)+$
                 ' '+lines[idx[j]].name
         endif else begin
            co='ff'x
            li=1
            factor=.1
            name=lines[idx[j]].name
         endelse
         tmp=min(abs(wave-lines[idx[j]].wave),ii)
         oplot,[1,1]*lines[idx[j]].wave, [1,sp[ii]],lines=li,co=co
         xyouts,lines[idx[j]].wave,sp[ii]*factor,name,orient=90
      endfor
;      stop
   endif
endfor
stop ;k loop
endfor

stop

return
end
