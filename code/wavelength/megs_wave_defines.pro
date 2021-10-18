;
; Define wavelength scale minima for MEGS channels
;
;min_megsb_L1_wave = 34.5d0  ;nm
;max_megsb_L1_wave = 105.5d0 ;nm
min_megsb_L1_wave = 33.0d0  ;nm
max_megsb_L1_wave = 107.0d0 ;nm

min_megsa1_L1_wave = 3.0d0 ;nm
max_megsa1_L1_wave = 39.0d0 ;nm (primary 5-19)

min_megsa2_L1_wave = 3.0d0 ;nm
max_megsa2_L1_wave = 39.0d0 ;nm (primary 16-38)

;megsa_1_cutoff_freq=6.4 ;nm ; version 1, 2
megsa_1_cutoff_freq=5.8D ;nm ; version 3
megsa_1to2_cutoff_freq = 17.24D	; nm

; June 26, 2014 bdt
; In level 2 processing this parameter is redefined
;	for the MA anomaly if needed by calling  set_megs_ab_cutoff_wave.pro
megs_AtoB_cutoff_freq = 37.0D	; nm

;
; Calculate number of bins for each channel
; for 0.01 nm sampling in L1 (or 100 bins/nm)
;
;num_megsa1_L1_spectral_elements = $
;	LONG((max_megsa1_L1_wave - min_megsa1_L1_wave + 1) * 100) ; 5 - 19 nm
;num_megsa2_L1_spectral_elements = $
;	LONG((max_megsa2_L1_wave - min_megsa2_L1_wave + 1) * 100) ;16 - 38 nm
;; all megs_a
;num_megsa_L1_spectral_elements = $
;	LONG((max_megsa2_L1_wave - min_megsa1_L1_wave + 1) * 100) ;5 - 38 nm

num_megsa1_L1_spectral_elements = $
	LONG((max_megsa1_L1_wave - min_megsa1_L1_wave) * 100) ; 5 - 19 nm
num_megsa2_L1_spectral_elements = $
	LONG((max_megsa2_L1_wave - min_megsa2_L1_wave) * 100) ;16 - 38 nm
; all megs_a
num_megsa_L1_spectral_elements = $
	LONG((max_megsa2_L1_wave - min_megsa1_L1_wave) * 100) ;5 - 38 nm
	
; megs_b
; Calculate number of bins for each channel
; for 0.02 nm sampling in L1 (or 50 bins/nm)
num_megsb_L1_spectral_elements  = $
	LONG((max_megsb_L1_wave  - min_megsb_L1_wave) * 50)  ;33.0 - 107.0 nm

num_megs_L2_spectral_elements = $
	LONG((max_megsb_L1_wave - min_megsa1_L1_wave ) * 50) ; 3 - 107 nm
	
;
; Define L1 MEGS wavelength scales
;
;megsa1_L1_wave = min_megsa1_L1_wave + 0.01*dindgen(num_megsa1_L1_spectral_elements)
;megsa2_L1_wave = min_megsa2_L1_wave + 0.01*dindgen(num_megsa2_L1_spectral_elements)
;megsa_L1_wave = range(min_megsa1_L1_wave,max_megsa2_L1_wave,invdelta=100,num_megsa_L1_spectral_elements)

megsa_L1_wave = min_megsa1_L1_wave + (0.01d0)*dindgen(num_megsa_L1_spectral_elements)
megsa1_l1_wave=megsa_l1_wave
megsa2_l1_wave=megsa_l1_wave

megsb_L1_wave  = min_megsb_L1_wave  + (0.02d0)*dindgen(num_megsb_L1_spectral_elements)
; range is limited to 8-digits of precision
;megsb_L1_wave  = range(min_megsb_L1_wave,max_megsb_L1_wave,invdelta=50,num_megsb_L1_spectral_elements)

;
; Define L2 wavelengths scale
;
;megs_L2_wave=    range(min_megsa1_L1_wave,max_megsb_L1_wave,invdelta=50,num_megs_L2_spectral_elements)
megs_L2_wave = min_megsa1_L1_wave + (0.02d0)*dindgen(num_megs_L2_spectral_elements) ; 3-106.98 inclusive

megs_L2_atob_cutoff_bin=MAX(WHERE(megs_l2_wave LT megs_atob_cutoff_freq)); Index of the last MEGS-A bin

