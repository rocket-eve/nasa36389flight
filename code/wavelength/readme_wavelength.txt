To get the wavelength image run

IDL> make_megsa_spectrum

This does a lot of stuff including comparisons with EVE flight data for
similar solar activity early in the mission, so some hand-tuning is
always necessary.
To check if all the code is being found use
IDL> .run make_megsa_spectrum
IDL> resolve_all
That should show most missing files. 

Known files that are needed include
do_ma_cm_calc1
fix_vc_offset
get_eve_ma_avg
make_megsa_spectrum
make_megsa_wave_img
make_sp_wave_cal
make_sp_wave_cal2
megsa_image_to_spectrum
megsb_defines
remove_vc_offset
sine_fun
sine2_fun
regrid_spectra
calc_2d_cm
calc_megs_fwhm
extract_megs_lines (seems broken)

Additional needed files
megsa1_solar_lines.dat (slit 1)
megsa_solar_lines.dat (slit 2)
reference_ma_spectrum_nrleuv_max.dat
regrid_spectra_mac.so (shouldn't be needed)
regrid_spectra_mac64.so

MEGS-B files needed
make_megsb_spectrum (main)
get_eve_mb_avg
megsb_image_to_spectrum
megsb_tuned_mask
plot_mb_lines
make_megsb_wave_img
make_megsb_stripe
megs_wave_defines

megsb_default_mask.sav
