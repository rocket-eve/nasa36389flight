function sine2_fun, x, p

; sum of two sine functions

return, p[0] + p[1] * ( sin( p[2] * x ) + sin( p[3] * x ) )
end
