function calc_2d_cm, img

iimg=double(reform(img))
cm=dblarr(2)
sum=total(iimg)
sx=total(iimg,2)
sy=total(iimg,1)
cm[0] = total( sx * dindgen(n_elements(iimg[*,0])) ) / sum
cm[1] = total( sy * dindgen(n_elements(iimg[0,*])) ) / sum
;stop

return,cm
end
