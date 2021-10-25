# nasa36353

This is the place to analyze flight data for NASA 36.353 SDO EVE
suborbital rocket measurements to produce one representative spectrum.
Much of it is based on previous analysis from other flights, SDO-EVE
and TIMED-SEE processing software. Most code is written in IDL, but the
atmospheric absorption code (not included) and the code to change
images into spectra (not included, uses call_external) are in fortran.
Intended for use by LASP personnel.

The general flow is to get each of these working.
Usually I setup the IDL_PATH in the rocket/number directory like this.
setenv IDL_PATH "+./:"${IDL_PATH}
That way I get the idl_lib and the <IDL_DEFAULT> paths. 
This needs read_netcdf.pro, and the EVE fits reader.

The first two convert raw images from Tom's sav files into calibrated spectra.
There are lots of steps in each.

1)
s=read_ma_36353()

2)
s=read_mb_36353()

3)
On evesci1, run the make_rkt..._atm_cor.pro procedure to capture MSIS model results.

4)
compare_ma_mb_36353
This one has a corresponding compare_ma_eve... and compare_mb_eve... that need
processed EVE data.
This merges all of the pieces together into one spectrum and writes the 
spectrum to a dat file.

