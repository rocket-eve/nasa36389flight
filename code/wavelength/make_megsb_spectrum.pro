function get_mb_rep_image, imglist_in

; need an even number of images to do this

; trim off the last image if there are an odd number of images
if n_elements(imglist_in[0,0,*]) mod 2 eq 1 then imglist_in=imglist_in[*,*,0:n_elements(imglist_in[0,0,*])-2]

;imglist = float(fix_vc_offset(imglist_in))
;; the shift by 4 pixels is no longer necessary here
imglist = imglist_in

; this is a dark-like correction, not needed
;for i=0L,n_elements(imglist[0,0,*])-1 do begin
;   imglist[*,*,i] = remove_vc_offset(reform(imglist[*,*,i]))
;endfor

img = median(imglist[*,*,0:*:2] < imglist[*,*,1:*:2],dim=3,/even)

; remove approx dark (looks pretty darn good already)
;stop

return,img
end

pro make_megsb_spectrum

!order=1
window,0,xs=1024,ys=512
loadct,39
device,decomp=0
if !p.color ne 0 then begin
   !p.background=!p.color
   !p.color=0
endif
co=independent_color()

@config36353

linefile='megsb_solar_lines.dat'

window,xs=1024,ys=512

;goto, do_img_restore

;restore,'./filtered/mbc.sav'
;mb=total(mbc,3)
;save,file='mbtot.sav',mb
;stop
;restore,'mbtot.sav' ;mb
;mb=mb/29.

; if the mb_corrected_imgs.sav file is not found, then
; we have to create it from the megs_b_raw_data.sav file
; in
; /eve_analysis/testing/analysis/Rocket_flight_anaylsis/rocket_#####/savesets/

;;;
;;; This creates the mb_corrected_imgs.sav file
;;;
;restore,'megs_b_raw_data.sav' ; megs_b structure
;restore,'36318_Flight_TM2_0_585_raw_bmegs.SAV' ; megs_b structure
;restore,'TM2_flight_36.318_raw_bmegs.SAV' ; megs_b structure
;restore,'data/flight_TM2_0_600_image_bmegs.sav' ; bmegs

workingdir = file_dirname(routine_filepath()) ; in wavelength dir
cd,workingdir

datapath = '../../data/'
restore,datapath+'rkt36'+numberstr+'_megsb_dark_particles_Oct_1_2021.sav'

megs_b = temporary(mb_no_spikes)

;stop
; ff at 8,9,10,11,12, 36, 80,81,82,83
;predarklist = [0,1,2,3,4,5,6,7, 13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,30,31,32, 41,42]
predarklist=[1,2,3,4,5]
; first solar image is 43
;sunlist=[45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75]
;sunlist=[13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37]
sunlist=[30,31,32,33,34,35,36,37,38,39,40]

; get a representative image
predark = get_mb_rep_image(megs_b[predarklist].image)
sun     = get_mb_rep_image(megs_b[sunlist].image)

images = (sun - predark ) ;> 1e-2
; good enough for finding a wavelength scale
;stop,'ready to save images to the file'

save,file='mb_corrected_imgs.sav',/compress, images ; contains scalar images variable

;;;
;;;
;;;
do_img_restore:

;restore,'best_mb.sav'
restore,'mb_corrected_imgs.sav' ;images
mb=images
goto,done_with_restore

;; read the data files
;;filelist=file_search('l0b_mb/MB*',count=count)
;;filelist=file_search('/Volumes/l0-b1/2010/099/MB__L0B_3_2010099_00*',count=count)
;stride=3
;while count mod stride ne 0 do count--
;filelist=filelist[0:count-1]
;tmp=fltarr(stride,2048,1024)
;imgs = fltarr(count/stride,2048,1024)
;for i=0L,count-1,stride do begin
;   for j=0L,stride-1 do begin
;      d=read_whole_fits(filelist[i+j < (n_elements(filelist)-1)],s)
;      ;tvscl,hist_equal(congrid(d.data.img,1024,512))
;      ;print,filelist[i+j]
;      ;stop
;      tmp[j,*,*]=d.data.img
;   endfor
;   med = median(tmp,dim=1)
;   tvscl,hist_equal(congrid(med,1024,512))
;   imgs[i/stride,*,*]=med
;endfor
;medimgs=imgs
;
;;darkfiles=file_search('/Volumes/l0-b1/2010/085/MB__L0B_1_2010085_200*',count=count)
;darkfiles=file_search('/Volumes/l0-b1/2010/094/MB__L0B_1_2010094_120*',count=count)
;stride=3
;while count mod stride ne 0 do count--
;darkfiles=darkfiles[0:count-1]
;tmp=fltarr(stride,2048,1024)
;dark = fltarr(count/stride,2048,1024)
;ccd_temp = fltarr(count)
;for i=0L,count-1,stride do begin
;   for j=0,stride-1 do begin
;      d=read_whole_fits(darkfiles[i+j < (n_elements(darkfiles)-1)],s)
;      ;tvscl,hist_equal(congrid(d.data.img < 2000,1024,512))
;      ;print,d.hdr.readout_mode
;      ;print,darkfiles[i]
;      ;stop
;      tmp[j,*,*]=d.data.img
;      ccd_temp[i+j] = d.hdr.ccd_temp
;   endfor
;   med=median(tmp < 2000,dim=1)
;   dark[i/stride,*,*]=med
;   ;tvscl,hist_equal(congrid(med,1024,512))
;endfor
;darkimg = median(dark,dim=1)

;; now remove darkimg from all data images
;for i=0L,n_elements(imgs[*,0,0])-1 do begin
;   imgs[i,*,*] -= darkimg
;   ; fudge
;   imgs[i,*,0:511] -= median(imgs[i,*,0:511])
;   imgs[i,*,512:1023] -= median(imgs[i,*,512:1023])
;   tvscl,hist_equal(congrid(reform(imgs[i,*,*])>0,1024,512))
;   ;stop
;endfor
;mb = (total(imgs,1) - 14) > 1e-10
;
;save,file='best_mb.sav',mb, imgs, darkimg, dark, ccd_temp, medimgs
;stop

done_with_restore:
tvscl,hist_equal(congrid(mb,1024,512))

stop

restore,'megsb_default_mask.sav'
zmask=default_mask
;zmask[0:1024,*]=shift(zmask[0:1024,*],0,10)

mb_orig=mb ;keep original just in case

megsb_tuned_mask, mb, zmask, zzmask

;add 5 pixels to each side, just for safety
 ; for 36.318, trim it down
for i=0L,2047 do begin
   gd=where(zzmask[i,*] eq 1,n_gd)
   if n_gd gt 1 then begin
      mmin=min(gd,max=mmax)
      ;zzmask[i, (mmin[0]-5) : (mmax[0]+5)] = 1b ; both sides
      zzmask[i, (mmin[0]+1) : (mmax[0]-1)] = 1b ;trim oen more off each side
   endif
endfor
zmask=zzmask
;stop

; DO NO DO THIS FOR 36.353
;zmask = shift(zmask,[0,26]) ;this adjusts the EVE mask to the rocket position
;; add one more pixel for 36.318
;zmask = shift(zmask,[0,6]) ;this adjusts the EVE mask to the rocket position

; for 36.318 the EVE mask is not close enough to match the rocket
; perhaps this has to do with the explosion from 36.300
; adjust the mask for 36.318
;  one side of the mask seems to be close, but the other is pretty far off

stripe = make_megsb_stripe(mb,zmask)
tvscl,congrid(stripe>0<10,1024,512)

;mb=(mb*zmask)>(1e-9)
mb=(mb*zmask)>(-.5)

;make_megsb_wave_img, mb, wimg, spwave, spout, file=linefile
;megsb_image_to_spectrum, mb, spwave, spflux,imgwave=waveimg,linefile=linefile,$
;                         imgmask=zmask
restore,'rkt36336_megsb_full_wave.sav' ; need an initial guess wavelength map
guess_imgwave = temporary(waveimg) ; 2048,1024
; tweak mask to match this other one that looks better
zmask = zmask and (guess_imgwave gt 1.)

; get the stripe image
stripeimg = make_megsb_stripe( mb, zmask )
thisimgwave = make_megsb_stripe( guess_imgwave*zmask, zmask ) ; 2048, 167
make_megsb_wave_img, stripeimg, thisimgwave, spwave, spflux, $
   file=linefile
waveimg = make_megsb_stripe(/undo, thisimgwave, zmask) ; convert back up to 2048x1024


;print,'saving rkt36336_mb_countspectrum.sav'
;save,file='rkt36336_mb_countspectrum.sav',/compress,spwave,spflux

plot,spwave,spflux>1,xs=1,ys=1,/ylog,xtit='Wavelength (nm)',ytit='counts/sec/slit height'

print,'saving rocket reference wavelengths'
description=systime(0,/utc)+' Fit from Raw TLM image file on macL4131 /Users/dlwoodra/idl/rocket/36'+numberstr+'/ with megsb_image_to_spectrum.pro'
save,file='rkt36'+numberstr+'_megsb_full_wave.sav',description,waveimg,/compress
stop
stop
stop


linewave = [36.8070, 49.9407, 58.4335, 62.9732, 79.0199, 94.9745, 97.2538, $
            97.7020, 102.5724, 103.1914]

get_eve_mb_avg, evesp

;; Michaels (Rachel's) spectrum
;;restore,'rocket_36290_combined_spectrum.sav' ; megs_irradiance, megs_wavelength

window,xs=1024,ys=800
!p.multi=[0,1,2]
!p.charsize=1.5
fitsp = spflux ; fill with correlation model
colors=[0,'fe'x,50]
step=5

for i=35,105,step do begin
   plot,spwave,spflux,xr=[i-2,i+step],xs=1,ys=1,ps=10,/ylog,yr=[1,1e5], $
        xtit='Wavelength (nm)', $
        ytit='Arbitrary (counts from rocket)', $
        tit='36.336 2018169 hour 19 UT (EVE is red)'
   legend_dlw,['rocket counts', 'EVE (max scale)', 'EVE ('+strtrim(step,2)+'nm linear corr)'],$
              /top,/left,textcolors=colors
   for j=0,n_elements(linewave)-1 do oplot,linewave[j]*[1,1],[1e-10,1e10],lines=1
   x=where(spwave gt i and spwave lt (i+step)<106.)
   scale = total(spflux[x]) / total(evesp[x])
   oplot,spwave,evesp*scale,co=colors[1],ps=10
   ; calculate a linear fit to model EVE to match the rocket
   r=poly_fit(evesp[x],spflux[x],1)
   print,i,r[0],r[1]
   fitsp[x] = evesp[x]*r[1] + r[0]
   oplot,spwave,fitsp, co=colors[2],ps=10
   plot,evesp[x],spflux[x],xtit='EVE L1 cps',ytit='rocket counts',ps=4,xs=1,ys=1
   oplot,evesp[x],fitsp[x],co=colors[2]
   stop
endfor

plot,spwave,spflux,/ylog,xr=[34,106],xs=1,yr=[1e-1,1e6],ps=10
oplot,spwave,fitsp,co=colors[2],ps=10
legend_dlw,['Rocket Counts','EVE linear fit counts'],textcolor=colors[[0,2]]

stop
stop
stop

goto,skip_h_compare

sumer=read_png('spectra_snaps/sumer920-940.png')
window,xs=n_elements(sumer[0,*,0]),ys=n_elements(sumer[0,0,*])
!order=0
tv,sumer,/true
!order=1

if !p.color ne 0 then begin
   !p.background=!p.color
   !p.color=0
endif
plot,spwave*10.,spflux/8.,/ylog,yr=[10,1e7],xr=[920,940],/noerase, $
     xmargin=[15.2,2.5],ymargin=[4.7,2.3],ys=5,xs=5,co=40
restore,'ch_sp_qs_0-1100nm.sav'
oplot,data.lambda,data.spectrum*(2e-8),co='f0'x,ps=10

stop

skip_h_compare:
plot,spwave,spflux,/ylog,yr=[10,1e7],xr=[32,109],ys=1,xs=1

lines=read_dat(linefile)

x=where(spflux gt 100) & lo=min(x) & hi=max(x)
lin=calc_megs_fwhm(spwave[lo:hi],spflux[lo:hi])

peaks=find_peaks(spwave,spflux)
stop
gidx=where(lin.width lt 0.5,n_gidx)
;for i=0,n_gidx-1 do print,lin[gidx[i]].w,lin[gidx[i]].width

gpeakrec={w:0.d,width:0.d,max:0.d,cm:0.d,chlinew:0.d,chion:''}
gpeaks=replicate(gpeakrec,n_elements(peaks))
gpeaks.w=peaks.w
gpeaks.max=peaks.max
gpeaks.cm=peaks.cm
gpeaks.width = calculate_gausswidth(spwave,spflux,peaks.w)

; now find closest matching chianti line
for i=0,n_elements(gpeaks)-1 do begin
   tmp=min(abs(gpeaks[i].w - lines.wave),idx)
   idx=idx[0] ;force to be a scalar
   gpeaks[i].chlinew = lines[idx].wave
   gpeaks[i].chion   = lines[idx].name
endfor


glinrec={w:0.d,pixel:0.d,max:0.d,width:0.d, chlinew:0.d0, chion:''}
glin=replicate(glinrec,n_elements(lin[gidx]))
glin.w     = lin[gidx].w
glin.pixel = lin[gidx].pixel
glin.max   = lin[gidx].max
glin.width = lin[gidx].width
; now find closest matching chianti line
for i=0,n_elements(glin)-1 do begin
   tmp=min(abs(glin[i].w - lines.wave),idx)
   idx=idx[0] ;force to be a scalar
   glin[i].chlinew = lines[idx].wave
   glin[i].chion   = lines[idx].name
endfor
save,file='mb_sp.sav',spwave,spflux,lines,lin,glin,gpeaks

for i=0,n_gidx-1 do $
   print,glin[i].w,glin[i].width,glin[i].chlinew,' ',glin[i].chion

; plot sections and label the lines that were found in glin
for i=30,105,5 do begin
   xrange=[i,i+6]
   plot,spwave,spflux,xr=xrange,xs=1,/ylog,ps=-4,yr=[10,1e7]
   for j=0L,n_elements(gpeaks)-1 do begin
      if gpeaks[j].w gt xrange[0] and gpeaks[j].w lt xrange[1] then begin
         str=strtrim(string(gpeaks[j].chlinew,form='(f8.3)'),2)+ $
             ' '+gpeaks[j].chion
         y=gpeaks[j].max
         oplot,[1,1]*gpeaks[j].w,[1,y],lines=1
         oplot,[1,1]*gpeaks[j].cm,[1,y],lines=2,co=co.green
         oplot,[1,1]*gpeaks[j].chlinew,[1,y],lines=2,co=co.blue
         xyouts,gpeaks[j].chlinew,y*1.5,orient=90,str
         xyouts,gpeaks[j].chlinew,y*0.1,orient=90,strtrim(string(gpeaks[j].width,form='(f8.3)'),2)
      endif
   endfor
;   for j=0L,n_elements(glin)-1 do begin
;      if glin[j].w gt xrange[0] and glin[j].w lt xrange[1] then begin
;         str=strtrim(string(glin[j].chlinew,form='(f8.3)'),2)+' '+glin[j].chion
;         y=glin[j].max
;         oplot,[1,1]*glin[j].w,[1,y],lines=1
;         oplot,[1,1]*glin[j].chlinew,[1,y],lines=2,co=co.blue
;         xyouts,glin[j].chlinew,y*1.5,orient=90,str
;         xyouts,glin[j].chlinew,y*0.1,orient=90,strtrim(string(glin[j].width,form='(f8.3)'),2)
;      endif
;   endfor
   stop
endfor
stop

plot_mb_lines,spwave,spflux,lines
stop

return
end
