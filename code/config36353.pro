; include file for 36.353

; to use, specify as @config36353 in the code
; in IDL these lines are then executed in place in that code the same as if it was typed into there
; include files ("batch" files in IDL) cannot support loops or any multiline commands or code blocks

; define constants used all over the place

version = '1_0'
numberstr = '353' ; 36.353
theyd = 2021252 ; sept 9, 2021 17:25:00 UT

;MEGS-A indices
dark1idx_a=[81,82,84,89,90] ; remove_megsa_spikes
ffidx_a=[133,134]
solaridx_a=[94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129]
dark2idx_a=[130,131,136,137]

;MEGS-B indices
dark1idx_b=[81,83,84,86,89,90] ; remove_megsb_spikes
ffidx_b=[133,134]
solaridx_b=[94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129]
dark2idx_b=[130,138,139]