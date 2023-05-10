pro get_eve_mb_avg, sp

; 36.258
;d = eve_read_whole_fits('MB__L1_2010123_18_004_01.fit.gz')

; 36.275
;d = eve_read_whole_fits('MB__L1_2011082_17_004_01.fit.gz')

; 36.286
;d = eve_read_whole_fits('MA__L1_2012175_19_004_01.fit.gz')

; 36.318
;d = eve_read_whole_fits('MB__L1_2016153_19_005_01.fit.gz')

; 36.336
;d = eve_read_whole_fits('MB__L1_2018169_19_006_02.fit.gz')

; 36.353
;d = eve_read_whole_fits('MB__L1_2021252_17_007_01.fit.gz')

; 36.389
d = eve_read_whole_fits('MB__L1_2023123_18_007_01.fit.gz')

dd = d.hdr001 ; data

; irradiance
;sp1 = total(d.hdr001.slit1_irradiance,2,/double) / double(n_elements(d.hdr001))
;sp2 = total(d.hdr001.slit2_irradiance,2,/double) / double(n_elements(d.hdr001))

; counts
;sp1 = total(d.hdr001.slit1_cps_per_pixel,2,/double) / double(n_elements(d.hdr001))
;sp2 = total(d.hdr001.slit2_cps_per_pixel,2,/double) / double(n_elements(d.hdr001))
sp = total(d.hdr001.cps_per_pixel,2,/double) / double(n_elements(d.hdr001))

;stop

return
end
