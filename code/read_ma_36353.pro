;+
; This function removes the dark from all of the MEGS-A images in the
; sequence. It uses linear interpolation from a pre-solar dark and
; post-solar dark reference image pair. It was tuned for the 36.258
; rocket flight in May 3, 2010 from WSMR.
; 
; :Params:
;    amegs: in, required, type=array of structures
;      This is the 
; :Returns:
;    A new array of structures is returned that has float type for the
;    image (with dark removed) and a new dark float image,
;    the total_counts tag is also updated to reflect real signal
;-
function remove_megsa_dark_36353, amegs, output

; determine a pre-observation dark and a post-observation dark
; goal is to linearly interpolate dark for each solar image

; change data type of image to float
imgdim=size(amegs[0].image,/dim)
newrec={time:0.d, pixel_error:0L, fid_index:0L, total_counts:0.d, $
  image:fltarr(imgdim), dark:fltarr(imgdim)}
output = replicate(newrec,n_elements(amegs)) ; create output structure

; copy data that is not changing
output.time = amegs.time
output.pixel_error = amegs.pixel_error
output.fid_index = amegs.fid_index

; DEFINE A DARK REGION MASK WHERE SUN LIGHT IS NOT DETECTED ON THE
; SHORT WAVELENGTH SIDE
darkmask=bytarr(2048,1024) + 1b ; all good
darkmask[0:3,*]=0b & darkmask[2044:2047,*]=0b ; remove edges
; what part should be used?, mask out everything else
; need a vertical strip on each half
darkmask[100:*,0:511]=0b ; keep 4-120 short wavelength side
darkmask[120:*,512:1023]=0b ; keep 4-120 short wavelength side
;stop
;stop
;stop







; define good image indices for dark images 
; good images do not show noise stripes across the width
; of the CCD
;preidx=[1,2,3,4,7,8,9]
;postidx=[49,54,55,57]
;preidx=[0,5,7,9,11,13,15,17,18,24,27,29]
;postidx=[68,69, 74,75,76]
preidx=[81,82,84,89,90]
postidx=[130,131,136,137]

; determine a pre-observation dark and a post-observation dark
; goal is to linearly interpolate dark for each solar image
predark = median(float(amegs[preidx].image),dim=3)
postdark= median(float(amegs[postidx].image),dim=3)

; prevent spikes in postdark
; plot postdark-predark and look for "normal" scatter
; the individual pixels should not be changing much
lobad=where(postdark-predark lt -6,n_lobad)
if n_lobad gt 0 then predark[lobad]=postdark[lobad]

; prevent spikes in predark
hibad=where(postdark-predark gt 5,n_hibad)
if n_hibad gt 0 then postdark[hibad]=predark[hibad]

; define a reference index (time) for the predark and postdark
; use time weighting to estimate custom dark for each solar image
pretime=preidx[-1] ; mean(preidx) ; 36.258
posttime=postidx[0] ;mean(postidx)
for i=0,n_elements(output)-1 do begin
  if i lt pretime then output[i].dark = predark else begin
     if i gt posttime then output[i].dark = postdark else begin
        ; normal case i > pretime and i < posttime
        ; interpolation is just a special combination of a weighted sum
        prefrac = (i-pretime) / (posttime-pretime)
        postfrac = 1. - prefrac
        output[i].dark = (predark*prefrac) + (postdark*postfrac)
     endelse
  endelse
endfor

; adjust row-by-row
for i=0,n_elements(output)-1 do begin
;;;;;  output[i].dark = postdark ; 9/8/20 This line was overriding the
;;;;;  weighting from the previous loop
  ; adjust dark to match row by row
  ; only allow up to +40 DN above dark and -20 DN below dark
  deltaimg = (darkmask*(amegs[i].image-output[i].dark))<20L > (-20L)
  for row=0L,1023 do begin
     ; use darkmask pixel differences, these are always dark
     ; except for flatfields
     diffarr = median(deltaimg[where(darkmask[*,row] eq 1),row],19)
     ; instead could we try a histogram-approach and select the peak?
     ; that would give the most likely, not the center value (median)
     ; probably needs more samples to determine accurately
     ; sticking with median
     gd=where(diffarr gt -20 and diffarr lt 20,n_gd)
     diff = mean(diffarr[gd]) ; exclude clipped portions
     ; skip flatfields
     if abs(diff) gt 19.9 then diff=0. ; just use normal dark
;if abs(diff) gt 1 then stop
     output[i].dark[*,row] += diff
  endfor
endfor

;stop
; subtract the dark everywhere
output.image = float(amegs.image) - output.dark
output.total_counts = total(total(output.image,1,/double),1,/double)

return,0
end

;+
;-
pro filterimgspike, arr, output_in

; needs at least 6 images for heavy filtering

output = output_in[arr]       ; only use specified subset 9/8/20 for 36.258 2010 rocket
  
imgdim=size(output[0].image,/dim)

; loop over only the specified subset of images in output
; calculate mean of each image
imgmean = mean(mean(output.image,dim=1),dim=1)

enable_heavy_filter = 0 ; only heavily filter dark
; dark means are around 0.04
; solar means are around 7.0
; flatfield means are over 6000.0

if mean(imgmean) lt 2 and n_elements(output) gt 5 then begin
   ; this is dark corrected dark data
   enable_heavy_filter=1
endif


if enable_heavy_filter eq 0 then begin
for i=0L,n_elements(output)-1 do begin

   next = (i+1) mod (n_elements(output)) ; wrap to compare last to first
   ; looks like each tap has a different noise threshold
   ; 0:511 (slit 2 with sam)
   bad=where(output[i].image gt output[next].image+sqrt(output[next].image>1.)*10.,n_bad)
   bot=output[i].image
   bot[bad]=output[next].image[bad]
   botspike=fltarr(imgdim)
   botspike[bad] = output[i].image[bad]

   ; 512:1023 (slit 1)
   bad=where(output[i].image gt output[next].image+sqrt(output[next].image>1.)*18.,n_bad)
   top=output[i].image
   top[bad]=output[next].image[bad]
   topspike=fltarr(imgdim)
   topspike[bad] = output[i].image[bad]
   
   output[i].image[*,0:511] = bot[*,0:511]
   output[i].spikes[*,0:511] = botspike[*,0:511]
   output[i].image[*,512:*] = top[*,512:*]
   output[i].spikes[*,512:*] = topspike[*,512:*]

endfor
endif
;stop

   
; dark heavy filter added for 36.258 2010 rocket
; this cannot filter consecutive strikes at the same pixel very well
tmpout = output
if enable_heavy_filter ne 0 then begin
   for i=0L,n_elements(output)-1 do begin
      next = (i+1) mod (n_elements(output))           ; wrap to compare last to first
      prev = (i-1)                                     ; minus 1 is OK since IDL wraps automatically
      tmp = (output[next].image + output[prev].image)*0.5 ; mean
      ;bad=where(output[i].image gt tmp+5., n_bad) ; clip to only 5 DN above mean to allow some normal noise
      ;if n_bad gt 0 then begin
      ;   tmpout[i].spikes[bad] = output[i].image[bad] ; copy bad pixels to spikes image
      ;   tmpout[i].image[bad] = tmp[bad]              ; replace bad pixels
      ;endif

      ; superior filtering uses a running 5-image median, temporal mean, and vertial 3 pixel median
      ; this removes waves in the dark images
      ; it might remove solar signal in solar images
      lo=(i-2) > 0
      hi=(lo+5) < (n_elements(output)-1L)
      lo = hi-5 ; re-adjust lo value if hi would be beyond the max
      ; can't do this next line with solar images
      ref = median(output[lo:hi].image,dim=3) < (tmp+5.) < (median(output[i].image,3)+5.) ; combine spatial median
      
      bad=where(abs(output[i].image - ref) gt 5., n_bad) ; clip to only 5 DN above mean to allow some normal noise
      if n_bad gt 0 then begin
         tmpout[i].spikes[bad] = output[i].image[bad] ; copy bad pixels to spikes image
         tmpout[i].image[bad] = ref[bad]               ; replace bad pixels
      endif

   endfor
   
   output = tmpout
endif

;; plot results
;for i=0,n_elements(output)-1 do begin
;   ; correct input index is arr[i] for input
;   plot,yr=[-20,20], output_in[arr[i]].image,tit='filterimgspike #'+strtrim(i,2)
;   oplot,output[i].image,co='fe'x
;; for solar, look at images
; tvscl,hist_equal(congrid(output_in[arr[i]].image,1024,512))
; stop
; tvscl,hist_equal(congrid(output[i].image,1024,512))
;   ;stop
;endfor

output_in[arr] = output ; copy back to input argument
return
end


;+
;-
function remove_megsa_spikes_36353, input, zzmask=zzmask

imgdim=size(input[0].image,/dim)
newrec={time:0.d, $
  image:fltarr(imgdim), spikes:fltarr(imgdim)}
output = replicate(newrec,n_elements(input)) ; create output structure

output.time = input.time
output.image = input.image

; first cut is to remove jumps up pairwise, comparing to next image
; if current image is larger than next image by sqrt(DN) then replace with next
; only work on "good" solar data

; first cut is to remove jumps up pairwise, comparing to next image
; if current image is larger than next image by sqrt(DN) then replace with next
; only work on "good" solar data

;identify types
;dark1idx=[0,3,5,7,9,11,13,15,16,18,22,24,26,27,29]
;ffidx=[2,25,70]
;solaridx=[31,32,33, 35,36, 38,39,40,41,42,43, 45,46,47,48,49,50,51,52,53,54, 56,57,58,59,60,61,62,63,64,65,66]
;dark2idx=[67,68,69, 71,72,73,74,75,76] ; only 68 has no waves

dark1idx=[81,82,84,89,90]
ffidx=[133,134]
solaridx=[94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129]
dark2idx=[130,131,136,137]

filterimgspike, dark1idx, output
filterimgspike, ffidx, output
filterimgspike, solaridx, output
filterimgspike, dark2idx, output

return,output
end

;+
; Use this function to prevent divide by zero errors for any value
; that is within tolerance of zero.
; Values less than tolerance are set to tolerance.
;
;-
function make_nozero_img, img, tolerance=tolerance

  if size(tolerance,/type) eq 0 then tolerance=1.d-12
  outimg=img
  bad=where(abs(outimg) lt tolerance,n_bad)
  if n_bad gt 1 then outimg[bad]=tolerance

return,outimg
end

;+
;-
pro make_ma_spectra36353, input, spectra, zmask1, zmask2, sens1img=sens1img, sens2img=sens2img, sens1_sp=sens1_sp, sens2_sp=sens2_sp, sens1errimg=sens1errimg, sens2errimg=sens2errimg, sens1err_sp=sens1err_sp, sens2err_sp=sens2err_sp, wimg1=wimg1, wimg2=wimg2

restore,'rkt36353_megsa_full_wave.sav'
; get l1 wavelength scale


wimg1=waveimg[*,512:*] ; slit 1
wimg2=waveimg[*,0:511] ; slit 2

; create a mean image to get initial values
; discover the solar line max

;ma=mean(input[20:30].image,dim=3) ; 36.336/290
tmp = reform(input.image[1628,256]) ; a 30.4 signal at each time
sm = smooth(median(tmp,3),7,/edge_trunc) ; a time smoothed value of 30.4
tmp = max(sm,maxidx) ; peak in time

ma=mean(input[maxidx-5:maxidx+5].image,dim=3) ; 36.258
;zero out sam area so it doesn't bias the results
zmask=byte(ma)*0b
zmask++

;; remove SAM
;zmask[300:600,0:511]=0

 ;center strip
;zmask[*,480:580]=0
zmask[*,440:515]=0

; zero out the slit 1  edge to help remove particle noise
zmask[*,860:*]=0
; zero out the slit 2  edge to help remove particle noise
zmask[*,0:70]=0
zmask[0:3,0:511]=0
zmask[2043:2047,512:1023]=0

; 36.258 SPECIAL MASK SHIFT
; MEGS-A IS SHIFTED BY ABOUT 70 PIXELS UPWARD RELATIVE TO 36.336
zmask = shift(zmask,0,70) ; 36.258

zzmask=zmask ; for returning
zmask1=zmask[*,512:*]
zmask2=zmask[*,0:511]

ma *= zmask
ma1=ma[*,512:*]
ma2=ma[*,0:511]


; make a reference spectrum
nozeroimg1 = make_nozero_img(ma1*zmask1,tol=1e-12)
;nozeroimg1 = ma1 * zmask1
;bad=where(abs(nozeroimg1) lt 1e-12,n_bad)
;nozeroimg1[bad]=1e-12

nozeroimg2 = make_nozero_img(ma2*zmask2,tol=1e-12)
;nozeroimg2 = ma2 * zmask2
;bad=where(abs(nozeroimg2) lt 1e-12,n_bad)
;nozeroimg2[bad]=1e-12

megsa_image_to_spectrum, nozeroimg2, sp2wave, sp2flux, imgmask=zmask2,imgwave=wimg2
megsa_image_to_spectrum, nozeroimg1, sp1wave, sp1flux, imgmask=zmask1,imgwave=wimg1

spectra_rec={time:0.d, w1:sp1wave, w2:sp2wave, sp1:sp1flux, sp2:sp2flux, sp1err:sp1flux*0., sp2err:sp2flux*0.}
spectra=replicate(spectra_rec, n_elements(input))
spectra.time=input.time

for i=0,n_elements(input)-1 do begin
   ;nozeroimg=reverse(input[i].image * zmask)
   nozeroimg = make_nozero_img(input[i].image*zmask,tol=1e-12)
   ;nozeroimg=input[i].image * zmask
   ;bad=where(abs(nozeroimg) lt 1e-12,n_bad)
   ;nozeroimg[bad]=1e-12

   img2 = nozeroimg[*,0:511]
   megsa_image_to_spectrum, img2, sp2wave, sp2flux, imgmask=zmask2,imgwave=wimg2, imgerr=sqrt(abs(img2)), sperr=sp2err

   img1 = nozeroimg[*,512:*]
   megsa_image_to_spectrum, img1, sp1wave, sp1flux, imgmask=zmask1,imgwave=wimg1, imgerr=sqrt(abs(img1)), sperr=sp1err

   spectra[i].w1 = sp1wave
   spectra[i].w2 = sp2wave
   spectra[i].sp1 = sp1flux
   spectra[i].sp2 = sp2flux
   spectra[i].sp1err=sqrt(total(sp1err^2,2,/double))
   spectra[i].sp2err=sqrt(total(sp2err^2,2,/double))

   plot,sp2wave,sp2flux,/ylog,yr=[.1,1e6],tit='Image #'+strtrim(i,2),ps=10
   oplot,sp1wave,sp1flux,co='fe'x,ps=10

   if i eq maxidx then begin
     p = image( hist_equal(congrid(input[i].image>0,1024,512)), $
      image_dimensions=[1024,512], dimension=[1024,512], rgb_table=3 )
     p.save,'megsa_'+strtrim(maxidx,2)+'_img.png',width=1024,height=512
     p.close
   endif

   heap_gc
endfor

if size(sens1img,/type) ne 0 then begin
   errimg=sens1errimg
   nozeroimg = abs(make_nozero_img(errimg,tolerance=1e-12))
   ;megsa_image_to_spectrum, sens1img, sp1wave, sens1_sp, imgmask=zmask1,imgwave=wimg1, imgerr=sqrt(sens1img), sperr=sens1_err_sp
   megsa_image_to_spectrum, sens1img*zmask1, sp1wave, sens1_sp, imgmask=zmask1,imgwave=wimg1, imgerr=nozeroimg, sperr=sens1err_sp
endif

if size(sens2img,/type) ne 0 then begin
   errimg=sens2errimg
   nozeroimg = abs(make_nozero_img(errimg,tolerance=1e-12))
   ;megsa_image_to_spectrum, sens2img, sp2wave, sens2_sp, imgmask=zmask2,imgwave=wimg2, imgerr=sqrt(sens2img), sperr=sens2_err_sp
   megsa_image_to_spectrum, sens2img*zmask2, sp2wave, sens2_sp, imgmask=zmask2,imgwave=wimg2, imgerr=nozeroimg, sperr=sens2err_sp
endif
stop
return
end

;+
; Need to tune for each rocket to try to correct bad images that are misassembled
;-
function fix_corrupted_image, amegs

; replace corrupted images in 36.353 (old way uses a median)
; (old way) amegs[95].image  = median(amegs[[94,96]].image,dim=3,/even)
tmp = amegs[95].image
; top and bottom rows are OK, but middle needs to be fixed
; have to fix top and bottom separately
rtmp = rotate(tmp,2) ; 180 deg rotation
; now need to shift to move VC to the right place
cn = 344
rtop = rtmp[*,0:511]
rtop = shift(rtop,cn) ; guess

rbot = rtmp[*,512:1023]
rbot = shift(rbot,-1*cn)
; examine then reassemble

new=tmp ; keep lowest 70 and highest 70 rows
rn = 70
new[*,rn:511] = rtop[*,rn:511]
new[*,512:1023-rn] = rbot[*,0:511-rn]
new[*,rn] = amegs[94].image[*,rn] ; copy from previous image
new[*,1023-rn] = amegs[94].image[*,1023-rn]

amegs[95].image = new

return, amegs
end
 




;+
; Custom rocket flight analysis for MEGS-A only
; 
; No parameters or keywords are used.
;
; :Examples:
;   IDL> s = read_ma_36353()
;-
function read_ma_36353

@config36353

window,0,xs=10,ys=10
wdelete
!p.color=0 & !p.background='ffffff'x & !p.charsize=1.5

; this expects a /code/ dir to contain this file and a /data/ dir at the same level as /code/ containing Tom's saesets

; change to this directory
workingdir = file_dirname(routine_filepath())
cd,workingdir

tomsMASaveFile = workingdir+'/../data/TM2_36353_Flight_MEGS-A_adata.sav' ; contains adata
;tomsMBSaveFile = workingdir+'/../data/TM2_36353_Flight_MEGS-B_bdata.sav' ; contains bdata

;ADATA           STRUCT    = -> <Anonymous> Array[140]

;IDL> help,adata,/str
;** Structure <253c128>, 5 tags, length=4194328, data length=4194328, refs=1:
;   TIME            DOUBLE          -809.89897
;   PIXEL_ERROR     LONG                -2
;   FID_INDEX       LONG            147630
;   TOTAL_COUNTS    DOUBLE       2.8983247e+09
;   IMAGE           UINT      Array[2048, 1024]

if file_test(tomsMASaveFile) ne 1 then begin
  print,'ERROR: cannot locate Toms MA save file '+tomsMASaveFile
  stop
endif
;if file_test(tomsMBSaveFile) ne 1 then begin
;  print,'ERROR: cannot locate Toms MB save file '+tomsMBSaveFile
;  stop
;endif

restore,tomsMASaveFile
; need to create amegs array of structures to match 36.290 adn 36.336
rec = {time:0.d, pixel_error:0L, fid_index:0, image:fltarr(2048,1024)}
;tmp = replicate(rec,n_elements(images[0,0,*]))

tmp = replicate(rec,n_elements(adata))
tmp.time = adata.time
tmp.pixel_error = adata.pixel_error
tmp.fid_index = adata.fid_index
tmp.image = float(adata.image)

; need data to be called amegs, fix corrupted images, too
amegs = fix_corrupted_image(tmp)

; keep backwards for now, need to fix VC offsets before reversing
; to put in wavelength order do amegs.image = reverse(temporary(images),1)
;amegs.image = images ; no reversal yet


;restore,tomsMBSaveFile
;description = info              ; info only stored in mb.sav
;reltime = info[0,*] ; seconds since launch (negative before)
reltime = amegs.time


amegs.time = reform(reltime)
heap_gc
; restores structure amegs with 77 images

; look at each image, and the difference with the previous one
xsize=1920 & ysize=1024
;window,0,xs=xsize,ys=ysize
;window,1,xs=800,ys=400
;window,2,xs=800,ys=400
for i=0,n_elements(amegs)-1 do begin
   if amegs[i].time lt 2090 then continue
   wset,0
   tmpimg = float(amegs[i].image)
   tmpimg[*,0:511] -= median(tmpimg[0:3,0:511])
   tmpimg[*,512:1023] -= median(tmpimg[2043:2047,512:1023])
   ;tvscl,hist_equal(congrid(amegs[i].image,xsize,ysize))
   tvscl,hist_equal(congrid(tmpimg,xsize,ysize))
;   device,decomp=0
;   tvscl,hist_equal(congrid(amegs[i].image mod 256, xsize, ysize))
   wset,1
   previdx = (i + n_elements(amegs) - 1) mod n_elements(amegs) 
   deltaimg = float(amegs[i].image) - float(amegs[previdx].image)
   !p.multi=0
   plot,deltaimg,title='index #'+strtrim(i,2)+'-#'+strtrim(previdx,2)+' T='+strtrim(amegs[i].time,2),yr=median(deltaimg)+[-1.,1.]*6.*stddev(deltaimg),ytit='Difference (DN)
   wset,2
   !p.multi=[0,1,2]
   plot,total(deltaimg,1)/2048.,ystyle=1,xstyle=1, $
      title='horizontal mean difference of #'+strtrim(i,2)+'-#'+strtrim(previdx,2),xtit='column',ytitle='mean DN/column'
   plot,total(deltaimg[*,0:511],2)/512.,ystyle=1,xstyle=1,$
      title='vertical mean difference of #'+strtrim(i,2)+'-#'+strtrim(previdx,2),xtit='row',ytitle='mean DN/row'
;   device,decomp=1
   oplot,total(deltaimg[*,512:1023],2)/512.,co='fe'x
   if i eq 95 then stop
end
;stop

;36.258
; 0 dark has vertical spike streak and frozen curve shape from power-up residual - discard (~9 DN off from final image)
; 1 dark shows weak horizontal bright/dim alternating patterns
; 2, 3 dark - discard, center values are decreasing as time goes forward and the detector settles down
; 4 dark shows similar horizontal ripple noise pattern - discard, center values are decreasing
; 5-7 dark, center values are decreasing
; discard all previous images
; 8 dark shows some larger period ripple noise pattern - ok? center values seem to have stabilized
; 9,10 dark
; 11 dark shows horizontal ripple pattern - discard
; 12-14 dark
; 15 dark shows some horizontal ripple pattern - discard
; 16,17 dark
; 18 dark shows some horizontal ripples - discard
; 19,20 dark
; 21 dark shows some horizontal ripples - discard
; 22,23 dark
; 24 dark has both horizontal and diagonal ripple pattern - discard
; 25, 26 dark has diagonal ripple patterns in slit 1 side
; 27 dark shows both horizontal and diagonal ripple pattern - discard
; 28 dark
; 29 dark shows some horizontal ripples - discard
; 30 dark
; 31 dark has diagonal ripple patterns in slit 1 side
; 32 dark shows both horizontal and diagonal ripple pattern - discard
; 33,34 dark has diagonal ripple patterns in slit 1 side
; 35 is partial RC with partial darks - discard
; 36 full RC
; 37 full RC has horizontal and diagonal ripple pattern
; 38, 39, 40 full RC
; 41, 42 partial dark discard
; 43-45 test pattern
; 46 partial test pattern/dark
; 47 dark shows both horizontal and diagonal ripple pattern - discard
; 48 dark has diagonal ripple patterns in slit 1 side
; 49 dark shows both horizontal and diagonal ripple pattern - discard
; 50,51 dark has diagonal ripple patterns in slit 1 side
; 52 dark shows both horizontal and diagonal ripple pattern - discard
; 53 dark has diagonal ripple patterns in slit 1 side
; 54 dark shows both horizontal and diagonal ripple pattern - discard
; 55,56 dark has diagonal ripple patterns in slit 1 side
; 57 dark shows both horizontal and diagonal ripple pattern - discard
; 58 dark has diagonal ripple patterns in slit 1 side
; 59 FF partially on
; 60 FF full, 61 FF is dimmer, 62 FF is dimmer, 63 partially on FF
; 64 dark has diagonal ripple patterns in slit 1 side
; 65 dark shows both horizontal and diagonal ripple pattern - discard
; 66 dark has diagonal ripple patterns in slit 1 side
; 67 dark shows both horizontal and diagonal ripple pattern - discard
; 68,69 dark has diagonal ripple patterns in slit 1 side
; 70,71 dark has diagonal ripple patterns in slit 1 side
; 72 dark shows both horizontal and diagonal ripple pattern - discard
; 73,74 dark has diagonal ripple patterns in slit 1 side
; 75 dark shows both horizontal and diagonal ripple pattern - discard
; 76 dark has diagonal ripple patterns in slit 1 side
; 77 dark shows both horizontal only near center and diagonal ripple pattern - ok?
; 78,79 dark has diagonal ripple patterns in slit 1 side
; 80 dark shows both horizontal and diagonal ripple pattern - discard
; 81,82 dark has diagonal ripple patterns in slit 1 side
; 83 dark shows both horizontal and diagonal ripple pattern - discard
; 84 dark has diagonal ripple patterns in slit 1 side
; 85 dark has diagonal ripple patterns in slit 1 side and some corruption near center from data loss
; 86 dark shows both horizontal and diagonal ripple pattern - discard
; 87 FF partially on
; 88 dark shows both horizontal and diagonal ripple pattern - discard
; 89,90 dark has diagonal ripple patterns in slit 1 side
; 91 dark shows both horizontal and diagonal ripple pattern - discard
; 92 first light smear - acquiring sun - discard
; 93 nearly full spectrum - acquiring sun, not centered - discard
; 94 sun centered full spectrum - horizontal pattern and diagonal pattern on slit 1
; 95 data loss - corrupted - FIXED in fix_corrupted_image
; 96 sun centered full spectrum - diagonal pattern on slit 1
; 97 sun centered full spectrum - horizontal pattern and diagonal pattern on slit 1 (slit 2 horizontal is mostly below the slit 2 spectrum toward edge)
; 98 sun centered full spectrum - diagonal pattern on slit 1
; 99 sun centered full spectrum - horizontal pattern and diagonal pattern on slit 1 (slit 2 horizontal is mostly above the slit 2 spectrum toward center)
; 100,101,102,103 sun centered full spectrum - diagonal pattern on slit 1
; 104 sun centered full spectrum - diagonal pattern on slit 1 - two huge particle strikes diagonally in slit 1
; 105,106,107 sun centered full spectrum - diagonal pattern on slit 1
; 108 sun centered full spectrum - horizontal pattern and diagonal pattern on slit 1 (slit 2 horizontal is mostly below the slit 2 spectrum toward edge)
; 109,110 sun centered full spectrum - diagonal pattern on slit 1
; 111 sun centered full spectrum - horizontal pattern and diagonal pattern on slit 1 (slit 2 horizontal is mostly below the slit 2 spectrum toward edge)
; 112 sun centered full spectrum - diagonal pattern on slit 1
; 113 sun centered full spectrum - horizontal pattern and diagonal pattern on slit 1 (slit 2 horizontal is mostly above the slit 2 spectrum toward center)
; 114,115 sun centered full spectrum - diagonal pattern on slit 1
; 116 sun centered full spectrum - horizontal pattern and diagonal pattern on slit 1 (slit 2 horizontal is mostly above the slit 2 spectrum toward center)
; 117 sun centered full spectrum - diagonal pattern on slit 1 (dark band under slit 1)
; 118 sun centered full spectrum - diagonal pattern on slit 1
; 119 sun centered full spectrum - horizontal pattern and diagonal pattern on slit 1 (slit 2 horizontal is in the center of slit 2)
; 120,121 sun centered full spectrum - diagonal pattern on slit 1
; 122 sun centered full spectrum - horizontal pattern and diagonal pattern on slit 1 (slit 2 horizontal is mostly below the slit 2 spectrum toward edge)
; 123,124 sun centered full spectrum - diagonal pattern on slit 1
; 125 sun centered full spectrum - horizontal pattern and diagonal pattern on slit 1 (slit 2 horizontal is mostly below the slit 2 spectrum toward edge)
; 126,127 sun centered full spectrum - diagonal pattern on slit 1 noticable dimmer solar signal
; 128 sun centered full spectrum - horizontal pattern and diagonal pattern on slit 1 (slit 2 horizontal is mostly above the slit 2 spectrum toward center)
; 129 sun centered full spectrum - diagonal pattern on slit 1 nearly all solar signal is gone from long wavelenghts of slit 2
; 130,131 dark has diagonal ripple patterns in slit 1 side
; 132,133,134 FF (135 is different)
; 136,137 dark has diagonal ripple patterns in slit 1 side
; 138,139 dark shows both lower frequency horizontal and diagonal ripple pattern - discard

preidx=[81,82,84,89,90] ; remove_megsa_dark
postidx=[130,131,136,137]

dark1idx=[81,82,84,89,90] ; remove_megsa_spikes
ffidx=[133,134]
solaridx=[94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129]
dark2idx=[130,131,136,137]

orig=amegs
; adjust locations by 4 pixels to fix real pixel locations
; shift columns to align lines across center (only one way will work)
; directions are opposite because readout is from opposite sides
; hard coded for TLBR
print,'INFO: moving VC to align top and bottom data'
for i=0,n_elements(amegs)-1 do begin
   amegs[i].image[*,0:511] = shift(amegs[i].image[*,0:511],-4,0)
   amegs[i].image[*,512:1023] = shift(amegs[i].image[*,512:1023],4,0)
   ; filter wild points
   bad=where(amegs[i].image ge (2.d^14),n_bad) ; larger than max DN
  if n_bad gt 0 then amegs[i].image[bad]=0 ; fill with zero
endfor

;put in wavelength order
print,'INFO: reversing images to put in wavelength order'
for i=0,n_elements(amegs)-1 do amegs[i].image=reverse(amegs[i].image)


print,'INFO: removing dark from all images'
status=remove_megsa_dark_36353(amegs, nodarkimg)

; now do particle filtering
print,'INFO: removing spikes'
nospikes=remove_megsa_spikes_36353(nodarkimg)

; image 104 needs work on slit 1

; save these results for Phil and the others
print,'saving intermediate data for analysis by others'
tmp=strsplit(systime(),/extract)
ts=tmp[1]+'_'+tmp[2]+'_'+tmp[4]
file='../data/rkt36353_megsa_dark_particles_'+ts+'.sav'
rec={time:0., image:fltarr(2048,1024)}
gd=where(nospikes.time gt -60,n_gd)
ma_no_spikes = replicate(rec,n_gd)
ma_no_spikes.time = nospikes[gd].time
ma_no_spikes.image = nospikes[gd].image
stop,'*** you really do not want to save this if it is not necessary***'
save,file=file, ma_no_spikes, /compress
ma_no_spikes=0b
heap_gc
print,'saved'
stop

return,-1
end


function temp

; image #44 doesn't seem to be filtered for spikes very well
bad=where(nospikes[44].image - nospikes[43].image gt 250,n_bad)
nospikes[44].image[bad] = nospikes[43].image[bad] ; manual fix
stop

;restore,'data/SURF_Sep10_sensitivity_megs_a1_fw0.sav' ; sensitivity includes gain
restore,'data/sensitivity_MEGSA1_second_order_2013.sav' ; sensitivity includes gain
old1=sens_uncorrected ;sensitivity
; fill old1
for i=0,2047 do begin
   gd=where(old1[i,*] gt 0 and old1[i,*] lt 100,n_gd,comp=comp)
   if n_gd gt 0 and n_elements(comp) gt 2 then begin
      old1[i,comp] = mean(old1[i,gd]) ; replace missing values with mean
   endif
endfor
bad=where(old1 lt 1e-8 or old1 gt 1e-2,n_bad)
old1[bad] = 1e-2
old1_sensitivity_error = old1*0.1 ;10% sensitivity_error


; new code 8/19 can remove the reverse
print,'INFO: getting slit 1 sensitivity'
a1sens = get_36353_sensitivity( /megsa1 )
a1sensitivity = a1sens[*,512:*] ; cut off



;; reverse so it is in wavelength order
; TODO: fix uncertainties (new ones are not good)
;; fold in an extra 10% beyond the 2010 uncertainties
a1senserr=reverse(old1_sensitivity_error[*,512:1023] * 1.1) > 1e-11 


restore,'data/36336sensitivity_MEGSA2_second_order_2017.sav' ; sensitivity from Brian 7/17/19
; MEGS-A2 save file from 36336 contains these items
;COMMENTS        STRING    = Array[4]
;SENS_CORRECTED1 FLOAT     = Array[2048, 1024] ; 
;SENS_CORRECTED2 FLOAT     = Array[2048, 1024]
;SENS_CORRECTED3 FLOAT     = Array[2048, 1024]
;SENS_UNCORRECTED
;                FLOAT     = Array[2048, 1024]
;WAVE            FLOAT     = Array[3600]

;adam308/evenetapp/store2/rocket_analysis_code/new_data/results/2013 > ls -l sensitivity_MEGS*
;-rw-rw-r-- 1 evesdp evesdp 16778840 Feb 19 19:07 sensitivity_MEGSA1_fw0_285MeV_2013.sav
;-rw-rw-r-- 1 evesdp evesdp 25182248 Feb 20 21:36 sensitivity_MEGSA1_second_order_2013.sav
;-rw-rw-r-- 1 evesdp evesdp 16778840 Feb 20 22:52 sensitivity_MEGSA2_fw0_380MeV_2013.sav
;-rw-rw-r-- 1 evesdp evesdp  8389880 Feb 21 16:51 sensitivity_MEGSA2_second_order_v2_2013.sav
;-rw-rw-r-- 1 evesdp evesdp 16778792 Feb 21 20:21 sensitivity_MEGSB_fw0_380MeV_2013.sav

old2 = reverse(sens_corrected1) ;sensitivity
; fill old2
for i=0,2047 do begin
   gd=where(old2[i,*] gt 0 and old2[i,*]lt 100,n_gd,comp=comp)
   if n_gd gt 0 and n_elements(comp) gt 2 then begin
      old2[i,comp] = mean(old2[i,gd]) ; replace missing values with mean
   endif
endfor
bad=where(old2 lt 1e-8 or old2 gt 1e-2,n_bad)
old2[bad] = 1e-2

old2_sensitivity_error = old2 * 0.1 ; 10%
;; update 7/19/19
;plot,old2[*,256],/ylog,yr=[1e-8,1e-4]
;oplot,sens_uncorrected[*,256],co='fe'x ; 2017 no 2nd order correction
;;oplot,sens_corrected1[*,256],co='fe00'x ; 2017 w/ 2009 2nd order
;;oplot,sens_corrected2[*,256],co='fe0000'x ; 2017 w/ 2009 183MeV
;;oplot,sens_corrected3[*,256],co='aaaa00'x ; 2017 w/ 2009 285MeV


print,'INFO: getting slit 2 sensitivity'
new2 = get_36353_sensitivity(/megsa2)

; filter
;;bad=where(new2 lt 1e-8 or new2 gt 1e-2,n_bad)
;bad=where(new2 lt 1e-8 ,n_bad)
;if n_bad gt 0 then new2[bad]=.001 ; make bad places insensitive

stop

; TODO: fix uncertainties (new ones are not good)
sensitivity_error = (old2_sensitivity_error * 1.1) > 1e-11
; fold in an extra 10% beyond the 2010 uncertainties

a2sensitivity = new2[*,0:511]  > 1e-8 ; 2048x512
a2senserr = sensitivity_error[*,0:511] > 1e-11 ; 2048x512

nospikes.image *= 0.1d ; divide by integration time 10 seconds
calimg=nospikes

stop

print,'INFO: integrating to create cps spectra'
make_ma_spectra36353, nospikes, spectra_cps, a1mask, a2mask, $
  wimg1=wimg1, wimg2=wimg2, $
  sens1img=a1sensitivity, sens1_sp=sens1_sp, sens1err_sp=sens1err_sp, $
  sens2img=a2sensitivity, sens2_sp=sens2_sp, sens2err_sp=sens2err_sp, $
  sens1errimg=a1senserr, sens2errimg=a2senserr

print,'INFO: applying sensitivity maps to cps images'
for i=0,n_elements(calimg)-1 do begin
   calimg[i].image[*,512:*] = calimg[i].image[*,512:*] * a1sensitivity
   calimg[i].image[*,0:511] = calimg[i].image[*,0:511] * a2sensitivity
endfor
print,'INFO: integrating to create calibrated spectra'
make_ma_spectra36353, calimg, spectra_cal, a1mask, a2mask, $
  wimg1=wimg1, wimg2=wimg2

a1sens2048=total(a1sensitivity*a1mask,2)/total(a1mask,2,/double) ; 2048
a2sens2048=total(a2sensitivity*a2mask,2)/total(a2mask,2,/double) ; 2048

a1senserr2048=total(a1senserr*a1mask,2)/total(a1mask,2,/double) ; 2048
a2senserr2048=total(a2senserr*a2mask,2)/total(a2mask,2,/double) ; 2048
; remove zeros
;cnt=0L
;while a1sens2048[cnt] lt 1e-12 do cnt++
;a1sens2048[0:cnt-1]=a1sens2048[cnt]
;a1senserr2048[0:cnt-1]=a1senserr2048[cnt]
;cnt=2047L
;while a1sens2048[cnt] lt 1e-12 do cnt--
;a1sens2048[cnt:*]=a1sens2048[cnt]
;a1senserr2048[cnt:*]=a1senserr2048[cnt]
a1sens2048[2043:2047]=a1sens2048[2042]

;cnt=0L
;while a2sens2048[cnt] lt 1e-12 do cnt++
;a2sens2048[0:cnt-1]=a2sens2048[cnt]
;a2senserr2048[0:cnt-1]=a2senserr2048[cnt]
;cnt=2047L
;while a2sens2048[cnt] lt 1e-12 do cnt--
;a2sens2048[cnt:*]=a2sens2048[cnt]
;a2senserr2048[cnt:*]=a2senserr2048[cnt]
a2sens2048[0:3]=a2sens2048[4]

; APPLY THE 2017/2009 SURF RATIO AS A DEGRADATION TERM TO ADJUST THE
; SENSITIVITY - S = S/deg
;restore,'rkt36290_degradation.sav' ; ma1deg, ma2deg, mbdeg
;stop
;a1sens2048 /= reverse(ma1deg)
;a2sens2048 /= reverse(ma2deg)
;a1senserr2048 ?
;a2senserr2048 ?
;stop

; wavelength
a1w2048=total(wimg1*a1mask,2)/total(a1mask,2) ; 2048
a1sens3600=interpol(a1sens2048,a1w2048,spectra_cps[0].w1,/nan)
; somtehing is wrong here...
a1senserr3600=interpol(a1senserr2048,a1w2048,spectra_cps[0].w1,/nan)

a2w2048=total(wimg2*a2mask,2)/total(a2mask,2) ; 2048
a2sens3600=interpol(a2sens2048,a2w2048,spectra_cps[0].w2,/nan)
a2senserr3600=interpol(a2senserr2048,a2w2048,spectra_cps[0].w2,/nan)

spectra=spectra_cps
for i=0,n_elements(spectra)-1 do spectra[i].sp1 = (spectra[i].sp1>.1)*a1sens3600
for i=0,n_elements(spectra)-1 do spectra[i].sp2 = (spectra[i].sp2>.1)*a2sens3600
;; DUE TO STRANGE EFFECTS, REPLACE SPECTRA with SPECTRA_CAL 10/2/19
stop
spectra = spectra_cal

k=3380 ;36.8 ; 3050 is 33.5
k17=1411 ; 17.11 nm
;plot,ps=-4,spectra_cps.sp2[k]/mean(spectra_cps[25:30].sp2[k]),yr=[-0.05,1.2],tit='Wavelength '+strtrim(string(spectra[0].w2[k],form='(f6.3)'),2)+' nm',ys=1,xr=[7,57],xs=1

; rel signal vs image #
msp = mean(spectra_cps[40:50].sp2,dim=2) ; mean "best" spectrum
k17col = 'fe'x ; red
plot,ps=-4,spectra_cps.sp2[k]/msp[k],yr=[-0.05,1.2],ys=1,xr=[27,68],xs=1,xtit='Image #',ytit='Relative signal'
oplot,ps=-1,spectra_cps.sp2[k17]/msp[k17],co=k17col
oplot,!x.crange,[0,0],lines=1
oplot,!x.crange,[1,1],lines=1
xyouts,/data,mean(!x.crange),.1,strtrim(string(spectra_cps[0].w2[k],form='(f6.3)'),2)+' nm'
xyouts,/data,mean(!x.crange),.2,strtrim(string(spectra_cps[0].w2[k17],form='(f6.3)'),2)+' nm',co=k17col
;oplot,[27,27],!y.crange,lines=1
stop

; rel signal vs time
time=spectra_cps.time
plot,ps=-4,time,spectra_cps.sp2[k]/msp[k],yr=[-0.05,1.2],ys=1,xr=[80,490],xs=1,xtit='L+ Time (Seconds)',ytit='Relative signal'
oplot,ps=-1,time,spectra_cps.sp2[k17]/msp[k17],co=k17col
xyouts,/data,mean(!x.crange),.1,strtrim(string(spectra_cps[0].w2[k],form='(f6.3)'),2)+' nm'
xyouts,/data,mean(!x.crange),.2,strtrim(string(spectra_cps[0].w2[k17],form='(f6.3)'),2)+' nm',co=k17col
oplot,!x.crange,[0,0],lines=1
oplot,!x.crange,[1,1],lines=1
stop

;comment='4/1/19 Switched from 380MeV on both slits to 285 and 140 for slit 1 and 2 respectively for version 3_6'
;comment='4/1/19 Switched variable names in sensitivity files, version 3_7'
;comment='7/25/19 MA1 sensitivity from 2017 w/ 2nd order from 2010, MA2 sensitivity from 2017 w/ 2nd order from 183MeV in 2009, version 3_8'
comment='09/04/20 Using sensitivity maps from 36.'+numberstr
save,file='rocket36'+numberstr+'_megsa_irr.sav',/compress,spectra,spectra_cps,a1sens3600,a1senserr3600, a2sens3600,a2senserr3600, a1sensitivity, a2sensitivity, a1sens2048, a2sens2048, comment
print,'saved rocket36'+numberstr+'_megsa_irr.sav'
stop

; no flight corrections, at earth
width = 800 & height = 600
;solar = lindgen(21)+14L         ; indices of solar measurements (14-35)
solar = lindgen(11)+40L ; 40-50
a1 = mean(spectra[solar].sp1 < 1.,dim=2)
a2 = mean(spectra[solar].sp2 < 1.,dim=2)
w = spectra[25].w1 ; same for both
gd=where(w gt 5.22 and w lt 38.39,comp=bad)
a1[bad]=1e-8 & a2[bad]=1e-8
gds1=where(w lt 20.5)
p = plot(w[gds1], a1[gds1], /stairstep, /ylog, $
    dimensions=[width,height], $
    xtitle='Wavelength (nm)',ytitle='Irradiance (W/m^2/nm)',$
    title='36.'+numberstr+' MA Slit Comparisons at Earth, no atm cor')
gds2=where(w gt 14) ; exclude SAM area
p =  plot(w[gds2], a2[gds2], color='red', /stairstep, /overplot)
p = text(/norm,.14,.9,'Slit 1',co='black')
p = text(/norm,.2,.9,'Slit 2',co='red')

; try to fix slit 1 wavelengths longer than 19 nm to match slit 2
keyw=[17.10,17.45,18.04,18.22,18.83,19.35,19.51,20.21]
r=a2/a1
index=lonarr(n_elements(keyw))
for i=0,n_elements(keyw)-1 do begin
   tmp=min(abs(w-keyw[i]),idx)
   index[i]=idx
endfor
; need to normalize r to be 1 at 17.1
coef=poly_fit(w[index],r[index]/r[index[0]],2,yfit=yfit)
print,'coef=',coef

corr=fltarr(n_elements(w)) + 1. ;initialize to 1.0
corr[1410:*] = poly( w[1410:*], coef ) > .0

gds1_c = where(w gt 17.2 and w lt 21)
p =  plot(w[gds1_c], (a1*corr)[gds1_c], color='blue', /stairstep, /overplot)
p = text(/norm,.14,.8,'Slit 1 polynomial corrected for 17-20 nm',co='blue')

plot,w,a2/a1,xr=[16,19],ps=10,ytit='A2/A1',xtit='Wavelength (nm)'
k = gaussian_function(10,width=51) & k /= total(k)
oplot,w,convol(a2,k)/convol(a1,k),co='fe'x
;plot,w[index],r[index]/r[index[0]],xs=1,ys=1
;oplot,w[index],yfit,co='fe'x

; approx values for a possible multiple gaussian fit
; 20.55=0.01, 21.57=0.001, 23.0=0.001, 26=.0001, 33=.00005

stop
stop
stop
return,1 ; the rest is unproven testing code to try to fix slit 1 to remove higher order lines

; try to guess second order spectrum in slit 1
model_ma1_sec_order_spec, w, a1, secorderspec

stop

ca1 = mean(spectra_cps[solar].sp1<(1.e5),dim=2)
ca2 = mean(spectra_cps[solar].sp2<(1.e5),dim=2)
w = spectra[25].w1 ; same for both
p = plot(w, ca1<1e5, /stairstep, $
    dimensions=[width,height], $
    xtitle='Wavelength (nm)',ytitle='Counts/Sec',$
    title='36.'+numberstr+' MA Slit Comparisons at Earth, no atm cor')
p =  plot(w, ca2, color='red', /stairstep, /overplot)
p = text(/norm,.14,.9,'Slit 1')
p = text(/norm,.2,.9,'Slit 2',co='red')
stop
stop
stop
stop


return,1
end
