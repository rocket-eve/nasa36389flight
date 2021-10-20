;+
;
; fix_vs_offset performs a vector shift to place the VC columns
; in the correct location. The virtual columns are not real and 
; the hardware uses these to replace real pixels. Since these are
; the first 4 pixels to be "read" from each line then the last 4 
; real pixels are lost on board.
;
; The top row on the rocket is from the left side so the top half
; of the assembled image needs to be shifted 4 pixels to the left.
; The bottom needs to be shifted 4 pixels to the right. This is 
; more complicated for different tap combinations, but the default
; tap combination on the rocket is TLBR (mode 1) and the rocket
; always flies with that.
;
; This code is vectorized to shift all tops at the same time and 
; all bottoms at the same time.
;
; :Params:
;    imgarr: in, required, type=imgarray
;       All primitive number types are OK. The dimentions are important.
;       DImensions must be 2048,1024,N in that exact order.
;
; :Example:
;  IDL> shuffled = fix_vc_offset( amegs.image )
;  Here the IDL evaluation of the array of structures containing
; 2d arrays is evaluated as [2048,1024,N] the way this code is
; expecting it.
;
;-
function fix_vc_offset, imgarr

;imgarr is [2048,1024,n]

dims = size(imgarr,/dim)
if dims[0] ne 2048 or dims[1] ne 1024 then begin
  print,'ERROR: fix_vs_offset did not revieve an array of size [2048,1024,N]'
  stop
endif

newimg = imgarr

newimg[*,0:511,*]    = shift(imgarr[*,0:511,*],[-4,0,0])
newimg[*,512:1023,*] = shift(imgarr[*,512:1023,*],[4,0,0])

return,newimg
end
