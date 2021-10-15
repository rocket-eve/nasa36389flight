pro get_eve_ma_avg, sp1,sp2

; 36.258
;d = eve_read_whole_fits('MA__L1_2010123_18_004_01.fit.gz')

; 36.275
d = eve_read_whole_fits('MA__L1_2011082_17_004_01.fit.gz')

;; 36.290
;d = eve_read_whole_fits('MA__L1_2013294_19_004_01.fit.gz')

dd = d.hdr001 ; data

; irradiance
;sp1 = total(d.hdr001.slit1_irradiance,2,/double) / double(n_elements(d.hdr001))
;sp2 = total(d.hdr001.slit2_irradiance,2,/double) / double(n_elements(d.hdr001))

; counts
sp1 = total(d.hdr001.slit1_cps_per_pixel,2,/double) / double(n_elements(d.hdr001))
sp2 = total(d.hdr001.slit2_cps_per_pixel,2,/double) / double(n_elements(d.hdr001))

stop

return
end
