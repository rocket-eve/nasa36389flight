;docformat = 'rst'

;+
; :Author:
;    Don Woodraska
;
; :Copyright:
;    Copyright 2019 The Regents of the University of Colorado.
;    All rights reserved. This software was developed at the
;    University of Colorado's Laboratory for Atmospheric and
;    Space Physics.
;
; :Version:
;  $Id$
;
;-

;+
; This is a wrapper for int_tabulated.
;
; :Params:
;   x_in: in, required, type=array
;     The x-values for each value of function f_in. These do not have to be uniformly spaced.
;   f_in: in, requlred, type=array
;     The y-values of the function. Must match number of x values.
;   lobins: in, required, type=array
;     Array of low bin edges in units of x_in
;   hibins: in, required, type=array
;     Array of high bin edges in units of x_in
;   status: out, optional, type=int
;     Return status code follows EXIS /mgunit convention. 1=good, 0=bad
; 
; :Returns:
;   The resultsarray is an array of integrated values (areas) for each element specified by lobins and hibins.
;   The units are units of x_in*f_in.
;   
; For integrating lyman-alpha over the EXIS FM1 and FM2 bandpasses (shift of -0.07nm) the ratio (1.0004) indicates 0.04% change.
; At 117.5 and 133.5 it is 0.2%. At 140 it is -0.8%, so brighter features as expected have smaller fractional error.
;
;-
function int_tabulated_array, x_in, f_in, lobins, hibins, status

  good_status=1
  bad_status=0

  if n_elements(x) ne n_elements(f) then begin
    print,'ERROR: int_tabulated_array - different lengths for x and f'
    return,bad_status
  endif
  if n_elements(lobins) ne n_elements(hibins) then begin
    print,'ERROR: int_tabulated_array - different lengths for lobins and hibins'
    return,bad_status
  endif

  tmpx = x_in
  tmpf = f_in

;  tmin = min(x_in) < min(lobins)
;  tmax = max(x_in) > max(hibins)
;  ; prepend and append zeroes far past the ends
;  ; special code to deal with boundaries
;  ; to reduce bad extrapolation, we append finite values outside provided range
;  if tmin lt 0 then tmin = tmin*2.d else tmin=0.d ; double lower bound or zero
;  if tmax gt 0 then tmax = tmax*2.d else tmax = 0.d ; double upper bound or zero
;  tmpx = [tmin, x_in, tmax] ; prepend and append points beyond input range
;  tmpf = [ 0.d, f_in, 0.d] ; prepend and append zeroes

  ; now insert all boundaries into the array so every bin can include the endpoints
  newx = [tmpx, lobins, hibins] ; join so all boundaries are part of the function
  newx = double(newx[ uniq(newx, sort(newx)) ]) ; sort unique values into increasing order
  nonnan = where(tmpf eq tmpf) ; skip nans
  newf = interpol(double(tmpf[nonnan]),double(tmpx[nonnan]),newx) ; linear interpolation to ensure there are values at bin edges

  resultsarray = dblarr(n_elements(lobins))

  for bin=0L,n_elements(lobins)-1 do begin
    ; find samples between lobins and hibins INCLUSIVE, edges are infinitely sharp
    gd = where( newx ge lobins[bin] and newx le hibins[bin], n_gd )
    if n_gd lt 2 then begin
      print,'ERROR: int_tabulated_array - not enough points in the bin to integrate - check lobins and hibins'
      stop
      return,bad_status
    endif
    resultsarray[bin] = int_tabulated(newx[gd], newf[gd], /double) ; arrays are already sorted
  endfor

  status = good_status
  return, resultsarray
end
