function remove_vc_offset, img

out = float(img)

; assume TLBR (default)
med = mean(median(img[2044:2047,1:510],3))
out[*,0:511] -= med
;for i=0,511    do out[*,i] = img[*,i] - mean( img[2044:2047,i] )

med = mean(median(img[0:3,513:1022],3))
out[*,512:1023] -= med
;for i=512,1023 do out[*,i] = img[*,i] - mean( img[0:3,i] )
;stop
return,out
end
