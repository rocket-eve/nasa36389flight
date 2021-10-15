function get_rep_image, imglist_in, no_fix=no_fix

; need an even number of images to do this

imglist = imglist_in
if keyword_set(no_fix) eq 0 then begin
   imglist = float(fix_vc_offset(imglist_in))
   ; not sure if the shift by 4 pixels is necessary here

   for i=0L,n_elements(imglist[0,0,*])-1 do begin
      imglist[*,*,i] = remove_vc_offset(reform(imglist[*,*,i]))
   endfor
endif

; if there are less than 3 images, then retun min image
if n_elements(imglist[0,0,*]) lt 3 then return,reform(imglist[*,*,0]<imglist[*,*,1])

; trim off last image if there are not an even number of images
if ((n_elements(imglist[0,0,*]) mod 2) eq 1) then imglist=imglist[*,*,0:-2] 
img = median(imglist[*,*,0:*:2] < imglist[*,*,1:*:2],dim=3,/even)

return,img
end

pro make_megsa_spectrum

!order=1
window,0,xs=1024,ys=512
loadct,39
;device,decomp=0
if !p.color ne 0 then begin
   !p.background='ffffff'x
   !p.color=0
endif

rktnum = '353' ; as in 36.353

; SLIT 1
linefile1='megsa1_solar_lines.dat'
; SLIT 2
linefile2='megsa_solar_lines.dat'

nrlmax=read_dat('reference_ma_spectrum_nrleuv_max.dat')

window,xs=1024,ys=512
;if file_test('ma_corrected_imgs.sav') eq 1 then goto, do_img_restore

; comment these two lines out to make a new saveset
;restore,'best_ma.sav' ; restores ma, imgs

; if the ma_corrected_imgs.sav file is not found, then
; we have to create it from the megs_a_raw_data.sav file
; in /eve_analysis/testing/analysis/Rocket_flight_analysis/rocket_#####/savesets/

;;;
;;; This creates the ma_corrected_imgs.sav file
;;;
;restore,'megs_a_raw_data.sav' ; megs_a structure
;restore,'TM2_flight_36.318_raw_amegs.sav'
;restore,'data/flight_TM2_0_600_image_amegs.sav'
;data = temporary(amegs)

workingdir = file_dirname(routine_filepath()) ; in wavelength dir
cd,workingdir

datapath = '../../data/'
restore,datapath+'rkt36'+rktnum+'_megsa_dark_particles_Sep_30_2021.sav'
data = temporary(ma_no_spikes)

predarklist=[1,2,3,4,5]
sunlist=[30,31,32,33,34,35,36,37,38,39,40]
postdarklist=[61,62]
; 

; get a representative image
predark  = get_rep_image(data[predarklist].image)
sun      = get_rep_image(data[sunlist].image)
postdark = get_rep_image(data[postdarklist].image)

images = (sun - predark) > 1e-2
; good enough for finding a wavelength scale

stop, 'ready to save images to the file'

save,file='ma_corrected_imgs.sav',/compress, images ; contains scalar images variable
;;;
;;;
;;;

do_img_restore:
restore,'ma_corrected_imgs.sav' ; images
ma=images + 1e-7

; ma is a "best" image to fit
heap_gc
;ma=reverse(ma) > .1
ma = (ma) > .1 ; already in wavelength order

device,decomp=0
tvscl,hist_equal(congrid(ma>1,1024,512))
stop

;kernel=[1.d0,2,3,4,3,2,1] ; 7 bin kernel
;kernel /= total(kernel)
;ama = convol(images,kernel,/edge_t,/center)
;tvscl,hist_equal(congrid(ama,1024,512))

;stop

;********
;stop

;zero out sam area so it doesn't bias the results
zmask=ma*0
zmask++

; remove SAM
zmask[300:600,0:511]=0

 ;center strip
;zmask[*,480:580]=0
zmask[*,440:515]=0

; zero out the slit 1  edge to help remove particle noise
zmask[*,860:*]=0
; zero out the slit 2  edge to help remove particle noise
zmask[*,0:70]=0
zmask[0:3,0:511]=0
zmask[2043:2047,512:1023]=0

ma *= zmask

; avg img (original way)
make_megsa_wave_img, ma[*,0:511], wimg2, sp2wave, sp2out, file=linefile2
make_megsa_wave_img, ma[*,512:*], wimg1, sp1wave, sp1out, file=linefile1,/slit1
!p.multi=[0,1,2]
plot,sp1wave,sp1out,/ylog,yr=[100,1e5],tit='Slit 1 36.'+rktnum,xtit='Wavelength (nm)'
plot,sp2wave,sp2out,/ylog,yr=[100,1e6],tit='Slit 2 36.'+rktnum,xtit='Wavelength (nm)'
stop
;; avg img
;make_megsa_wave_img, ama[*,0:511], wimg2, sp2wave, asp2out, file=linefile2
;make_megsa_wave_img, ama[*,512:*], wimg1, sp1wave, asp1out, file=linefile1,/slit1
;stop

waveimg = ma
waveimg[*,0:511]=wimg2 ;the 30.4 side with sam
waveimg[*,512:*]=wimg1 ;the short wavelength side

device,decomp=1
plot,sp1wave,sp1out>10,/ylog,xr=[15,25],yr=[1,1e6]
oplot,sp2wave,sp2out>10,co='fe'x

; oplot some lines
;lw=[ 6.9632, 7.5034,  9.3923, 9.6121, 9.8116, 10.0576, 10.3566, 10.5208, 12.7666, 13.1240, 14.8402, 15.2154, 15.4162, 17.1073, 
lw=[ 6.9661, 9.3923, 9.6121, 10.0576, 10.3566, 10.5208, 12.7666, 13.1240, 14.8402, 15.2154, 15.4162, 17.1073, $
  17.7240, 17.4532, 18.0401, 18.8232, 25.63175, 28.415, 30.3783, 33.541, 36.076, 36.8076]
for i=0,n_elements(lw)-1 do oplot,lw[i]*[1,1],[1,1e10],lines=1

; overplot the NRLEUV spectrum scaled for comparison
oplot,nrlmax[0,*],nrlmax[1,*]*1e9,co='ff0000'x,ps=10


stop
plot,sp1wave,sp1out>10,/ylog,xr=[25,37],xs=1,yr=[1,1e5]
oplot,sp2wave,sp2out>10,co='fe'x
for i=0,n_elements(lw)-1 do oplot,lw[i]*[1,1],[1,1e10],lines=1
; overplot the NRLEUV spectrum scaled for comparison
oplot,nrlmax[0,*],nrlmax[1,*]*1e9,co='ff0000'x,ps=10
stop

get_eve_ma_avg, eve_sp1, eve_sp2
; scale to 17.1 area on slit 1
tmp17=where(sp1wave gt 17 and sp1wave lt 17.2)
sf1=total(sp1out[tmp17])/total(eve_sp1[tmp17])
sf2=total(sp1out[tmp17])/total(eve_sp2[tmp17])
sf0=total(sp1out[tmp17])/total(sp2out[tmp17])
tmpn17=where(nrlmax[0,*] gt 17 and nrlmax[0,*] lt 17.2)
sfn=total(sp1out[tmp17])/total(nrlmax[1,tmpn17])

save,file='rkt36'+rktnum+'_countspectrum.sav', sp1wave, sp2wave, sp1out, sp2out, lw, nrlmax, eve_sp1, eve_sp2
print,'saved countspectrum'

!p.charsize=1.5
step=2 ; nm
colors=[0UL,'fe'x,'aaaa'x,'f00000'x,'f000'x]
for j=5,35,step do begin
   plot,sp1wave,sp1out>10,/ylog,xr=[j,j+step],xs=1,yr=[1,1e5],ps=10
   oplot,sp2wave,sp2out*sf0>10,co=colors[1],ps=10
   for i=0,n_elements(lw)-1 do oplot,lw[i]*[1,1],[1,1e10]
   ; overplot the NRLEUV spectrum scaled for comparison
   oplot,nrlmax[0,*],nrlmax[1,*]*sfn,co=colors[2],ps=10

   x=where(sp2wave ge j and sp2wave lt j+step)
   r=poly_fit(eve_sp1[x],sp1out[x],1) ; linear
   oplot,sp2wave,eve_sp1*r[1] + r[0],co=colors[3],ps=10 ; blue

   r=poly_fit(eve_sp2[x],sp2out[x]*sf0,1) ; linear
   oplot,sp2wave,eve_sp2*r[1] + r[0],co=colors[4],ps=10 ; green

   ;oplot,sp2wave,eve_sp1*sf1,co=colors[3],ps=10 ; blue
   ;oplot,sp2wave,eve_sp2*sf2,co=colors[4],ps=10 ; green
   legend_dlw,['Slit1','Slit2','NRLEUV','EVE Slit1 2013294','EVE Slit2 2013294'],$
              textcolor=colors
   ;stop
endfor

do_ma_cm_calc1,ma

stop

; This code looks messed up

;; correct bad parts of both fits
;; short wavelengths for slit 2 and long wavelengths for slit 1
;; choose representative rows from each half
;
;;fix short wavelength side first
;waveimg2=waveimg
;; overwrite fitted wavelengths for the bad part of slit 2
;tmp=min(abs(waveimg[*,300]-waveimg[*,750]),crossover)
;; or just pick a crossover point
;;for i=0,crossover-1 do begin
;for i=0,700 do begin
;   waveimg2[i,0:511]=waveimg2[i,750]
;endfor
;;for i=crossover,2047 do begin
;for i=1100,2047 do begin
;   waveimg2[i,512:*]=waveimg2[i,300]
;endfor
;
;stop

;for tr=250,350 do begin
;   ;   675,775
;   ;pick corresponding row from the bottom
;   br=425+tr
;   ; use slit 1 (512:1023) to fix 0:680 on slit 2 (0:511)
;   med_offset=median(waveimg2[*,tr] - waveimg2[*,br])
;   waveimg2[0:680,tr] = waveimg2[0:680,br]+med_offset
;endfor
;
;;fix long wavelength side
;for br=675,775 do begin
;   tr=br-425
;   ;med_offset=median(waveimg2[680:*,tr] - waveimg2[680:*,br])
;   med_offset=(waveimg2[999,tr] - waveimg2[999,br])
;   waveimg2[1000:2047,br] = waveimg2[1000:2047,tr]-med_offset
;endfor
;
;;fill in the in-between spots
;; use vertical interpolation
;vert=dindgen(1024)
;;vg=[vert[251:349],vert[676:774]] ; vertical good (where spectra fall)
;vg=[vert[200:425],vert[625:850]] ; vertical good (where spectra fall)
;vb=[vert[0:199],vert[426:624],vert[851:1023]] ; vertical bad (between and outside)
;for i=0,2047 do begin
;;   xpos = where(waveimg2[i,*] eq waveimg[i,*],n_xpos,comp=comp,ncomp=ncomp)
;;   if n_xpos ne 0 and ncomp ne 0 then $
;;      waveimg2[i,xpos] = interpol(waveimg2[i,comp],vert[comp],vert[xpos])
;   waveimg2[i,vb] = interpol(waveimg2[i,vg],vert[vg],vert[vb])
;endfor
;;stop

;;recalculate sp2out
;wimg2 = waveimg2[*,0:511]
;make_megsa_wave_img,ma[*,0:511],wimg2, sp2wave,sp2out
;sp2out[0:100]=.1
;
;;recalculate sp1out
;wimg1 = waveimg2[*,512:1023]
;make_megsa_wave_img,ma[*,512:1023],wimg1, sp1wave,sp1out,/slit1
;sp1out[0:100]=.1
;
;stop
;waveimg = waveimg2
;
;; try to tweak wavelength based on CM of each line
;tweak_ma_cm_wave, ma, waveimg, waveimg2
;waveimg2 = waveimg
;stop

description=systime(0,/utc)+' Fit to each slit separately 36.'+rktnum+' data on macL4131 using /Users/dlwoodra/idl/rocket/36'+rktnum+' with make_megsa_spectrum.pro, note that slit 2 is row 0 to 511, and slit1 is row 512 to 1023.'
save,file='rkt36'+rktnum+'_megsa_full_wave.sav',description,waveimg,/compress
print,'---'
print,'saved rkt36'+rktnum+'_megsa_full_wave.sav'
print,'---'
   device,decomp=1
   plot,sp1wave,sp1out,xs=1,/ylog,yr=[100,1e8],/ys, ps=10, $
        tit='MA Rough Spectrum from image sum',xr=[5,20];xr=[30,39]
   oplot,sp2wave,10*sp2out,co='f0'x,ps=10
   oplot,nrlmax[0,*],10.*max(sp1out)*nrlmax[1,*]/max(nrlmax[1,*]),co=180,ps=10
   xyouts,6,1e7,'Slit 1'
   xyouts,7,1e7,'Slit 2',co='f0'x
   xyouts,6,1e6,'NRLEUV max',co='180'

   plot,sp1wave,sp1out,xs=1,/ylog,yr=[100,1e8],/ys, ps=10, $
        tit='MA Rough Spectrum from image sum',xr=[20,40];xr=[30,39]
   oplot,sp2wave,10*sp2out,co='f0'x,ps=10
   oplot,nrlmax[0,*],10.*max(sp1out)*nrlmax[1,*]/max(nrlmax[1,*]),co=180,ps=10

lines1=read_dat(linefile1)
lines2=read_dat(linefile2)
save,file='ma_sp.sav',sp1wave,sp1out,sp2wave,sp2out,lines1,lines2
stop
stop
stop

; compare to EVE spectra on the rocket day (36.318 on 2016153
; near 19 UT)
; NO MEGS-A data is available from EVE after the MEGS-A anomaly in 2014/146

;get_eve_ma_avg, eve_sp1, eve_sp2
; scale to 17.1 area on slit 1
tmp17=where(sp1wave gt 17 and sp1wave lt 17.2)
sf1=total(sp1out[tmp17])/total(eve_sp1[tmp17])
sf2=total(sp1out[tmp17])/total(eve_sp2[tmp17])
sf0=total(sp1out[tmp17])/total(sp2out[tmp17])
stride=2
for i=5,37,2 do begin
   plot,sp1wave,sp1out,xs=1,/ylog,yr=[1e1,1e5],ys=1,xr=[i,i+stride],ps=10
   oplot,sp2wave,sp2out*sf0,co=colors[1],ps=10 ; red
;   oplot,sp2wave,eve_sp1*sf1,co=colors[2],ps=10 ; blue
;   oplot,sp2wave,eve_sp2*sf2,co=colors[3],ps=10 ; green
   for j=0,n_elements(lw)-1 do oplot,lw[j]*[1,1],[1,1e10]
   legend_dlw,['Slit1','Slit2','EVE Slit1','EVE Slit2'],$
              textcolor=colors[0:3]
   if !p.multi[0] eq 0 then stop   
endfor

;print,'.c to plot spectra with lines'
;stop
;rl=[7.2312, 8.6914, 13.0941, 17.1073, 17.7240, 18.8232]
;loc=[206., 309., 605., 859.0, 897.0, 964.0]
;rl=[ 14.8402, 17.1073, 17.7240, 18.8232, 30.3783, 36.8076]
rl=[ 14.8402, 15.2154, 15.4162, 17.1073, 17.4532, 17.7240, 18.0401, 28.415,    30.3783,  36.076,  36.8076]
;loc=[206.,    309.,    605.,   859.0,   897.0,   964.0, 1620.75, 1957.]
;loc=[750.,    892.0,   912.0,   928.0,  1648,     1984.]
loc=[750.193,  774.014,  786.955,  891.368,  912.321, 928.852, 948.128, 1542.963, 1648.506, 1946.797, 1983.863];pre2010091
; day 2010/099
loc=[750.72589,774.50574,787.46119,891.65981,912.79418,929.31186,948.63751,1543.7008,1649.1427,1947.4931,1984.5597]

;rocket 36290
rl=[   14.8402,   15.2154,   15.4162,   17.1073,   17.4532,   17.7240,   18.0401, 19.5119, 23.733, 24.303,    28.415,   30.3783,  36.076,  36.8076] ; lines

;loc=[725.     , 749.     , 763.,      867.11467, 888.50421, 904.95026, 924.47446,   1014., 1262.,  1293,       1522, 1628.1273,    1927, 1965]
loc=[725.     , 749.     , 763.,      865.02275, 887.292,   903.187,   923.022,     1014,  1262,   1293,       1519.294, 1625.483, 1924., 1962.]
;rl=[17.107, 30.3783, 36.8076] ;slit 2
;loc=[862.25, 1623.75, 1960] ;slit 2
coef2=poly_fit(loc,rl,2,yfit=yfit2)
print,coef2 ;use coef1 as first guess at wavelength scale in megs_sp_wave_cal
stop
stop
stop

;tweak_ma_cm_wave,ma,waveimg,waveresult
;
;stop


sp=total(ma[*,512+210:512+220],2)
;w=megs_sp_wave_cal(sp,err,status,file=linefile1,/slit1,param=param)
w=megs_sp_wave_cal2(sp,err,status,file=linefile1,/slit1,param=param)
;plot_ma_lines,w,sp,lines1
sp2=total(ma[*,512-220:512-210],2)

rl=[6.9632, 7.5034, 9.4012, 9.6122, 10.3566, 14.8402, 17.1073, 17.4532, 17.7240, 18.0401]
;loc=[224.,  261.,     391.,  406.,    458.,    749.,   890.0,   912.0,   928.0,  948.0]
loc=[183.,  221.,     353.,  367.,    419.,    713.,   856.0,   877.0,   894.0,  913.0]
coef1=poly_fit(loc,rl,2,yfit=yfit1)
print,coef1 ;use coef1 as first guess at wavelength scale in megs_sp_wave_cal

stop
;plot_ma_lines,sp1wave,sp1out,lines1
;plot_ma_lines,sp2wave,sp2out,lines2


stop
sp1out = 200.*sp1out / max(sp1out) ;make a relative spectrum
;sp1out = 70.*sp1out / max(sp1out) ;make a relative spectrum
; try getting extract_megs_lines to work
sp1out[0:50]=sp1out[51]*.9
sp1out[3599]=sp1out[3598]*.9

x=where(sp1out gt 1.e-6) & lo=min(x) & hi=max(x)<(3595);hi=max(x)<(1809)
;linesfound = extract_megsa_lines( sp2out )
lines_f1=calc_megs_fwhm(sp1wave[lo:hi],sp1out[lo:hi])
gidx=where(lines_f1.width lt .175)
gausswidth = calculate_gausswidth(sp1wave[lo:hi],sp1out[lo:hi],lines_f1[gidx].w)
help,lines_f1
for i=0L,n_elements(gidx)-1 do begin
   print,lines_f1[gidx[i]].w,lines_f1[gidx[i]].width,' half',lines_f1[gidx[i]].w/2.
endfor
stop
stop
stop
x=where(sp2out gt 1.e-6) & lo=min(x) & hi=max(x)
lines_f2=calc_megs_fwhm(sp2wave[lo:hi],sp2out[lo:hi])
stop
stop
stop
;

;trial
sp=total(ma,2)
w=megs_sp_wave_cal(sp,err,s,file=linefile2)
plot,w,sp
stop

;make_megsa_wave_img, ma[*,0:511],wref, spwave, spflux,file=linefile2
megsa_image_to_spectrum, ma[*,0:511], sp2wave, sp2flux, line=linefile2

stop



stop

return
end
