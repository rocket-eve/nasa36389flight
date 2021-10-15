pro do_ma_cm_calc1, img

pmulti_orig=!p.multi

mac=(img-20.)>(1e-9) ;suppress noise less than 20 DN
; first do 30.4
cm304=fltarr(2)
cm171_2=cm304
cm171_1=cm304
cm368=cm304

;lw=[ 6.9632, 7.5034,  9.4012, 9.6121, 9.8116, 10.0576, 10.3566, 10.5208, 12.7666, 13.1240, 14.8402, 15.2154, 15.4162, 17.1073, $
; 17.4532, 17.7240, 18.0401, 18.8232, 30.3783, 36.8076]
lw=[ 14.8402, 15.2154, 15.4162, 17.1073, $
     17.4532, 17.7240, 18.0401, 28.415, 30.3783, 36.076, 36.8076]
; long wavelength side
;x=[     750L,     774,    787,      891, $
;         912,     929,    948,     1543,  1648,    1947, 1984]
x=[     725L,     749,    763,      867, $
         889,     905,    925,     1522,  1628,    1927, 1965]

halfwidth=3

;slit2 side (long wavelengths)
for i=0L,n_elements(lw)-1 do begin
   lo=x[i] - halfwidth
   hi=x[i] + halfwidth
   print,lw[i],'CM x & y =',calc_2d_cm(mac[lo:hi,0:511])+[lo,0]
   ;plot,lindgen(hi-lo+1) + lo,total(mac[lo:hi,0:511],2)
   ;stop
endfor
return

   cm304=calc_2d_cm(mac[1624:1634,0:511])
   cm368=calc_2d_cm(mac[1963:1967,0:511])
   cm171_2=calc_2d_cm(mac[863:871,0:511])
   cm171_1=calc_2d_cm(mac[863:871,512:*])

!p.multi=[0,1,2]
cm304[0,*] += 1624.
cm368[0,*] += 1963.
cm171_2[0,*] += 863.
cm171_1[0,*] += 863.
cm171_1[1,*] += 512.

plot,cm304[0,*],ps=-4,yr=[long(min(cm304[0,*])),ceil(max(cm304[0,*]))], $
     tit='Slit 2 30.4 nm CM-x',xtit='image #',ytit='pixel #'
plot,cm304[1,*],ps=-4,yr=[long(min(cm304[1,*])),ceil(max(cm304[1,*]))], $
     tit='Slit 2 30.4 nm CM-y',xtit='image #',ytit='pixel #'
stop

plot,cm368[0,*],ps=-4,yr=[long(min(cm368[0,*])),ceil(max(cm368[0,*]))], $
     tit='Slit 2 36.8 nm CM-x',xtit='image #',ytit='pixel #'
plot,cm368[1,*],ps=-4,yr=[long(min(cm368[1,*])),ceil(max(cm368[1,*]))], $
     tit='Slit 2 36.8 nm CM-y',xtit='image #',ytit='pixel #'
stop

plot,cm171_2[0,*],ps=-4,yr=[long(min(cm171_2[0,*])),ceil(max(cm171_2[0,*]))], $
     tit='Slit 2 17.1 nm CM-x',xtit='image #',ytit='pixel #'
plot,cm171_2[1,*],ps=-4,yr=[long(min(cm171_2[1,*])),ceil(max(cm171_2[1,*]))], $
     tit='Slit 2 17.1 nm CM-y',xtit='image #',ytit='pixel #'
stop

plot,cm171_1[0,*],ps=-4,yr=[long(min(cm171_1[0,*])),ceil(max(cm171_1[0,*]))], $
     tit='Slit 1 17.1 nm CM-x',xtit='image #',ytit='pixel #'
plot,cm171_1[1,*],ps=-4,yr=[long(min(cm171_1[1,*])),ceil(max(cm171_1[1,*]))], $
     tit='Slit 1 17.1 nm CM-y',xtit='image #',ytit='pixel #'
stop

plot,cm304[1,*]-mean(cm304[1,*]),ps=-4
oplot,cm171_2[1,*]-mean(cm171_2[1,*]),ps=-4
oplot,cm171_1[1,*]-mean(cm171_1[1,*]),ps=-4
stop
stop
!p.multi=pmulti_orig
return
end
