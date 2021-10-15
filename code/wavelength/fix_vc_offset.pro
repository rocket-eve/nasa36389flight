function fix_vc_offset, imgarr

;imgarr is [2048,1024,n]

newimg = imgarr

newimg[*,0:511,*]    = shift(imgarr[*,0:511,*],[-4,0,0])
newimg[*,512:1023,*] = shift(imgarr[*,512:1023,*],[4,0,0])

return,newimg
end
