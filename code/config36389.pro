; include file for 36.389

; to use, specify as @config36389 in the code
; in IDL these lines are then executed in place in that code the same as if it was typed into there
; include files ("batch" files in IDL) cannot support loops or any multiline commands or code blocks

; define constants used all over the place

version = '1_0'
numberstr = '389' ; 36.389
theyd = 2023123 ; may 3, 2023 18:30:00 UT
humandatestr='May 3, 2023'
au = 0.98419335 ; from lisird lasp_vsop87_1au_correction_PT1M
; call this interactively to get 1-AU factor from LISIRD
; s=get_lisird_data(dataset='lasp_vsop87_1au_correction_PT1M',mintime='2023-05-03T18:30:00',maxtime='2023-05-03T17:45:00',/jd)
; Earth is far sun, so need to increase irradiance
; The irradiace at 1au is the measured irradiance divided by au

apogeekmstr = '293.98' ; revisit
apogeesecstr = '280.05' ; check (279?)

; used in msis00e
ft7 = 164.3 ; adjusted
ft7a = 152.0                    ; TODO - revisit ft7a (81-day avg) when known
;more penticton_radio_flux.csv | tail -183 | awk -F',' '{print $2}' | awk -f ~/bin/avg.awk
fap =2.0 ; Potsdam - TODO - revisit (https://kp.gfz-potsdam.de/en/data)

sza=31.05 ; deg TODO: revisit

; 90 deg roll at T+370s
; solar observing T+150-470s

;MEGS-A indices
;dark1idx_a=[81,82,84,89,90] ; remove_megsa_spikes
;ffidx_a=[133,134]
;solaridx_a=[94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129]
;dark2idx_a=[130,131,136,137]

dark1idx_a = [5,13,14,15,16,17,20,21,22,23,24,27,28,29,30,31,34,35,36,37,38,42,43,44,45,46,49,50,51,52,57,58,61,62,65,66,67,68,69,70,74,75,76,77,78,82,83,84,90,91,92,93,94,95,99,100,101,102,103,107,110,111,115,116,117,118,119,123,124,125,126,127,128,132,133,134,135,136,140,141,142,143,144,147,148,149,150,151,155,156,157,158,159,174,178,179,180,181,182,193,194,195,196,200,201,202,203,204, 209,210,211,216,217,218]

ffidx_a=[187,188,189]

solaridx_a=[221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,$
            ; start of roll (sam change)
            247,248,249,$
            ; roll completed
            250,251,252,253,254,255,256] ; 257 has motion/dim
; for 90deg only use 250-256

dark2idx_a=[265,270,271,272,273,274,275,276,280,281,282,283,284,285,286,287]

;MEGS-B indices
;dark1idx_b=[81,83,84,86,89,90] ; remove_megsb_spikes
;ffidx_b=[133,134]
;solaridx_b=[94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129]
;dark2idx_b=[130,138,139]

;MEGS-B indices for 36.389
dark1idx_b=[2,3,6,7,8,10,11,14,15,16,17,18,21,22,23,24,25,29,30,31,32,33,36,37,38,39,40,43,44,45,46,47,48,52,54,55,61,62,63,64,68,69,70,71,72,76,77,78,84,87,88,89,93,94,95,96,101,102,103,104,105,110,111,117,118,119,122,126,127,128,129,130,134,135,136,137,138,142,143,144,145,146,150,151,152,153,154,158,159,160,161,173,174,175,176,181,182,183,184,195,196,197,198,202,203,204,205,206,210,211,212,213,215,216,217,218] ; needs to be divisible by 3?
; 212, 213, 216 are suspect

ffidx_b=[187,188,189,190,260,261]

dark2idx_b=[257,258,263,264,265,266,267,268,269,272,273,274,275,276,277,278,279,283,284,285,286,287]

solaridx_b=[221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,$
                                ; start of roll
           246,247,248,249,250,251,252,253,254,255]
; for 90deg only use 248-255

;59 corrupted/data loss
;79 change T=-1298, corruption
;85 change, corruption
;108 change T=-1008, corruption
;112 change T=-968, corruption
;113 corrupted/shifted
;120 corrupted/shifted upsidedown
;163 change to reverse clock partial T=-458
;167 RC off partial dark
;169 partial TP and dark
;172 TP off partial
;186 partial FF
;191 FF off
; launch at 209
; post launch dark jumps around from image to image
;214 partial FF
; 219 first partial solar spectrum moving
; 246 90deg roll start
; 256 T=471 last partial solar, shutter door close at T+475
; 259 partial ff
; 262 partial ff

; strong temperature dependent pixel at [1772,510] (and vertical)
; and at [1900,510] (and vertical)
; and [1591,510] (and vertical) [1592,510] (and vertical)
