pro regrid_spectra, wave_in, flux_in, n_spec, n_in, $
                    wg_in, em_in, ea_in, wave_out, flux_out, n_out, $
                    wg_out, em_out, ea_out, acor, dcor, acoro, dcoro, status
;+
; NAME:
;   regrid_spectra.pro
;
; PURPOSE:
;   Provide an interface to use the FORTRAN subroutine regrid_spectra.
;
; CATEGORY:
;   L2_EGS
;
; CALLING SEQUENCE:
;   regrid_spectra, wave_in,flux_in, n_spec, n_in, $
;                    wg_in, em_in, ea_in, wave_out, flux_out, n_out, $
;                    wg_out, em_out, ea_out, acor, dcor, acoro, dcoro, status
;
; INPUTS:
;     wave_in(n_in,n_spec): wavelengths for each spectrum
;     flux_in(n_in,n_spec): flux in irradiance units (W/m^2) for spectra
;     n_spec: number of spectra to regrid
;     n_in: number of points (wave,flux pairs) in each spectrum
;     wg_in(n_in,n_spec): byte array where flux data is good=1 or bad=0
;     em_in(n_in,n_spec): error (measurement precision of flux)
;     ea_in(n_in,n_spec): error (total accuracy of flux)
;     wave_out(n_spec): wavelengths for output
;     n_out: number of output wavelength bins
;     acor(n_in,n_spec): a 2d array to be regridded
;     dcor(n_in,n_spec): a 2d array to be regridded
;
; OPTIONAL INPUTS:
;     none
;
; KEYWORD PARAMETERS:
;     none
;
; OUTPUTS:
;     flux_out(n_out,n_spec): trapezoid integrated flux in wave_out
;     wg_out(n_out,n_spec): byte array where output data is good=1 or bad=0
;     em_out(n_out,n_spec): output bin error (measurement precision)
;     ea_out(n_out,n_spec): output bin error (total accuracy)
;     acoro(n_out,n_spec): regridded 2d array
;     dcoro(n_out,n_spec): regridded 2d array
;     status: 0 if OK, -1 if error occurred causing low flux in a bin
;
; OPTIONAL OUTPUTS:
;     none
;
; COMMON BLOCKS:
;     none
;
; SIDE EFFECTS:
;     speed increase over the original IDL routine average_egs_irradiance
;
; RESTRICTIONS:
;     output bins must be sorted in increasing order with no large gaps
;     wave_out must be higher resolution than wave_in
;      output bins cannot have more than one input bin (wave_out must
;      be smaller than wave_in)
;     each flux is considered to be independent of neighbor bins
;      (if data is correlated, underestimate of errors will occur)
;     IDL version 5.4 or newer is required to run this procedure
;      (because auto_glue is not supported until 5.4)
;     inputs and output are quietly forced to proper types for CALL_EXTERNAL
;
; PROCEDURE:
;     uses CALL_EXTERNAL to interface to FORTRAN subroutine
;
; MODIFICATION HISTORY:
;     10-5-01 Don Woodraska Original File Creations.
;
; $Log: regrid_spectra.pro,v $
; Revision 6.0  2002/09/12 15:31:29  dlwoodra
; update to 6.0
;
; Revision 5.1  2002/09/12 15:00:35  dlwoodra
; update from main
;
; Revision 5.20  2002/09/06 23:21:15  see_sw
; commit of version 5.0
;
; Revision 4.0  2002/05/29 18:08:45  see_sw
; Release of version 4.0
;
; Revision 3.0  2002/02/01 18:54:36  see_sw
; version_3.0_commit
;
; Revision 1.1  2001/10/31 22:13:25  dlwoodra
; initial commit
;
;
;idver='$Id: regrid_spectra.pro,v 6.0 2002/09/12 15:31:29 dlwoodra Exp $'
;-

if !version.release lt 5.4 then begin
    print,''
    print,' IDL 5.4 or higher is REQUIRED to execute this procedure.'
    print,' CALL_EXTERNAL uses /auto_glue keyword'
    print,''
    print,' FATAL ERROR ENCOUNTERED IN REGRID_SPECTRA.PRO'
    print,''
    stop,'STOPPED'
endif

dllpath = file_dirname(routine_filepath())+path_sep() ; same path as this code

case strlowcase(!version.os) of
   'linux' : begin
      dll = dllpath + 'regrid_spectra_lnx664.so'
      entry = 'regrid_spectra_'
   end
   'darwin' : begin
        ;dll = dllpath + 'regrid_spectra_mac.so'
        ;entry = 'regrid_spectra__' ;g77 appends extra underscore to image
        entry = 'regrid_spectra_' ;gfortran
        ; x86 no longer maintained for mac
        ;if strcmp(!version.arch,'x86_64') then begin
        ;   dll = dllpath + 'regrid_spectra_mac64.so'
        ;endif
        if strcmp(!version.arch,'arm64') then begin
           dll = dllpath + 'regrid_spectra_mac64.so'
        endif
    end
    'sunos' : begin
        ;
        ; Use 32-bit or 64-bit shareable object library?
        mbits   = ((!version.release ge 5.4) ? !version.memory_bits : 32)
        dll     = dllpath + string(mbits,"('regrid_spectra',i2,'.so')")
        entry   = 'regrid_spectra_'
    end
endcase

; quietly force inputs to proper type
wave_in  = float(TEMPORARY(wave_in))
flux_in  = float(TEMPORARY(flux_in))
n_spec   = long(n_spec)
n_in     = long(n_in)
wg_in    = byte(TEMPORARY(wg_in))
em_in    = float(TEMPORARY(em_in))
ea_in    = float(TEMPORARY(ea_in))
wave_out = float(TEMPORARY(wave_out))
n_out    = long(n_out)
acor     = float(TEMPORARY(acor))
dcor     = float(TEMPORARY(dcor))

;quietly force outputs to proper type
flux_out = TEMPORARY(fltarr(n_out,n_spec))
wg_out   = TEMPORARY(bytarr(n_out,n_spec))
em_out   = flux_out
ea_out   = flux_out
acoro    = flux_out
dcoro    = flux_out
status   = 0L


result=CALL_EXTERNAL(dll, entry, $
       wave_in, flux_in, n_spec, n_in, $
       wg_in, em_in, ea_in, wave_out, flux_out, n_out, $
       wg_out, em_out, ea_out, acor, dcor, acoro, dcoro, status,/auto_glue)

return
end
  
