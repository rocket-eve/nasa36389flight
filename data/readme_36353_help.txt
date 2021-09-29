Instructions for Data for NASA 36.353
----------------------------------------------------------
9/11/2021  Tom Woods

Flight Information
---------------------
Launched on Sept. 9, 2021 (DOY 252) at 17:25:00 UT  (11:25:00 MDT)
Solar observations are from T+150s to T+480s.
Apogee is at T+278s of 291 km.

B8 flare peaked at 17:16 UT and XRS-B level was down to B4 - B3 range during flight.
F10.7 = 102.

IDL SETUP
---------
The alias rocket_idl, rocket_data need to be defined before running IDL.
These definitions can live in your .cshrc startup file.

alias rocket_idl 'setenv IDL_PATH "+$home/Dropbox/minxss_dropbox/code/rocket_real_time:+$home/Dropbox/minxss_dropbox/code/dualsps:+$home/Dropbox/minxss_dropbox/code/minxss_library:<IDL_DEFAULT>"'

setenv rocket_data /Users/twoods/Dropbox/minxss_dropbox/rocket_eve_docs/36.353/TM_Data

TM Flight Data Files are stored in $rocket_data/Flight/TM1 and $rocket_data/Flight/TM2

Rocket TM Read Code lives mostly in
	/Users/twoods/Dropbox/minxss_dropbox/code/rocket_real_time/TM_Read_Code/

Terminal Window input to start up idl:
USER%  rocket_idl
USER%  idl

************************************************************************************************
**************************    IDL SAVE SET FILES FOR 36.353 FLIGHT ON SEPT 9, 2021 *************
**************************    Located in Dropbox/minxss_dropbox/rocket_eve/36.353/ *************
**************************    Solar observing between T+150s and T+470s            *************
************************************************************************************************

TM2 Files
---------
MEGS-A CCD images from T-810s to T+580s:  ./TM_Data/Flight/TM2/TM2_36353_Flight_MEGS-A_adata.sav
	Restores adata structure array.  Solar measurements are from T+150s to T+480s.

MEGS-B CCD images from T-810s to T+580s:  ./TM_Data/Flight/TM2/TM2_36353_Flight_MEGS-B_bdata.sav
	Restores bdata structure array.  Solar measurements are from T+150s to T+480s.

TM1 Files
---------

Analog Monitors:  ./TM_Data/Flight/TM1/analogs_36353_flight_adata.sav

ESP data:  ./TM_Data/Flight/TM1/esp_36353_flight_edata.sav
	Restores edata structure array.  Solar measurements are from T+150s to T+480s.

MEGS-P data:  ./TM_Data/Flight/TM1/megs_p_36353_flight_pdata.sav
	Restores pdata structure array.  Solar measurements are from T+150s to T+480s.
	This also includes MEGS-A and MEGS-B analog data.

DualSPS XRS+X55: ./TM_Data/Flight/TM1/dualsps_xrs_xdata_x55_messages.sav
	Restores packets: xdata, x55, messages
	Solar Data are for TIME_SINCE_ON for 950s to 1310s
	Dark Data are good for TIME_SINCE_ON for 100s to 800s
	SPS-1 visible quad: saturated - needs ND filter added or reduced op amp gain
	SPS-2 diodes: very low signal - needs increased op amp gain
	    sps2_data[0]=XRS-B,  sps2_data[1]=LY-Alpha,  sps2_data[2]=XRS-C, sps2_data[3]=Dark

OWLS EUVOP:  no useful science data (saturated diodes): ./TM_Data/Flight/TM1/dualsps_owls_odata.sav
	SPS-1 visible quad: needs ND filter added or reduced op amp gain
	SPS-2 diodes: need filters added


************************************************************************************************
**************************   PROCESSING / PLOTTING CODE IN IDL  ********************************
************************************************************************************************

The primary raw binary data files have full TM packets with a time stamp.
These WSMR-processed files are called LOG files.

Tom has processed these files (using instructions below) and extracted the data and saved into
*.dat binary files and also as *.sav IDL savesets for each instrument.
Others do not need to reprocess the Binary files.

NOTE:
IDL save set (*.sav) files are better to use than the binary data files due to new binary file
issues between different type computers in how they consider integers as either 32-bit or 64-bit integers.

As end user, see the following QUICKLOOK instructions for the data (*.dat) files.

**********   WARNING: while most code has been updated for 36.353 processing, there still might
be some code that works with the old rocket flights but not yet updated for 36.353 (2021 flight).

-----------------------------------------------------------
-----------------------------------------------------------
QUICKLOOK Instructions - How to plot / read the rocket data
-----------------------------------------------------------
-----------------------------------------------------------

MEGS CCD Images - movie & read single image
-------------------------------------------
IDL>  movie_raw_megs, MEGS_FILE_NAME, 'A' or 'B', info=info

This procedure is written to be able to read single image at a time:
IDL>  aimage = 8  ; index into images in the file
IDL>  movie_raw_megs, MEGS_FILE_NAME, 'A', image=aimage
;  on return, "aimage" is the 8th image

MEGS_FILE_NAME = 36353_Flight_TM2_TTT_TTT_raw_Xmegs.sav

where        X = "a" or "b"
and
       TTT-TTT = -200_100   for dark and flat field data - 29 images
               = Ê100_460   for primary solar data (near apogee) - 34 images
               = Ê460-585   for more dark and flat field data at end of flight - 12 images


ESP Data - plot & read data
---------------------------
IDL>  plot_esp, 'TM1_36353_Flight.log_esp.dat', esp_data
 or
IDL>  plot_esp, 'esp_36353_flight_edata.sav'

	Channel    Dark-1    Dark-2      Apogee   Visible   Signal (Light-Dark)
	ESP # 1     41.42     40.80       52.60      0.49       11.49
	ESP # 2     27.12     26.38      163.43     -0.44      136.68
	ESP # 3   1165.62   1189.25     1176.62     -8.56       -0.82  ; known bad channel
	ESP # 4     42.80     42.60      113.42      0.03       70.72
	ESP # 5     33.17     32.58       86.42     -0.15       53.54
	ESP # 6     50.60     50.58       97.28     -0.03       46.69
	ESP # 7     40.17     39.62      100.54     -0.19       60.64
	ESP # 8     44.10     43.78      647.78      0.21      603.85
	ESP # 9     33.20     32.72      339.47     -0.13      306.50

	QUAD X =  -1.62 arc-minutes   ; this is more about flare location than about sun center
	QUAD Y =  -9.77 arc-minutes

MEGS-P Data - plot & read data
---------------------------
IDL>  plot_megsp, 'TM1_36353_Flight.log_pmegs.dat', pmegs_data
 or
IDL>  plot_megsp, 'megs_p_36353_flight_pdata.sav'

	  Channel      Dark-1    Dark-2      Apogee  Signal (Light-Dark)
	MEGS-P # 1     27.17     27.05       32.27        5.15
	MEGS-P # 2     31.65     30.98       31.30       -0.01			; Dark Channel

MEGS-P packet also has the CCD analog monitors.  For plotting those in addition to the
photometer data, use the /all option.

IDL>  plot_megsp, 'TM1_36353_Flight.log_pmegs.dat', pmegs_data, /all


Dual-SPS Data - plot & read  data
---------------------------
IDL>  plot_dualsps, 'TM1_36353_Flight.log_sps_xrs.dat', 'SPS', data=data

Dual-SPS Read All Packets
-------------------------
IDL>  data = read_dualsps_file( 'TM1_36353_Flight.log_sps_xrs.dat', messages=messages, x55=x55, /verbose )
or
IDL>  restore, 'dualsps_xrs_data_x55_messages.sav'  ; restores DATA, X55, MESSAGES packets

Analog Monitors Data - plot & read data
----------------------------------------
IDL>  plot_analogs, 'TM1_36353_Flight.log_sps_analogs.dat', analog_data
 or
IDL>  plot_analogs, 'analogs_36353_flight_adata.sav'


---------------------------------------------------
---------------------------------------------------
1.  Instructions for Extracting Data from Binary Files (LOG files, full packets)

	***** Tom has done this for the FLIGHT data, so others do NOT have to do this.
---------------------------------------------------
---------------------------------------------------

LOG Binary File Names
-----------------
*.LOG are TM files with full packets  with time stamps.

TM1_36353_Flight.log  is TM#1 file (3 Mbps) for Analog monitors, ESP, XPS, MEGS-P, Dual-SPS units

TM2_36353_Flight_ver2.log is TM#2 file (10 Mbps) for CCD data (MEGS-A & B)

NOTE: Be sure to use the /cd option for LOG files for read_tm1.pro and read_tm2.pro.

TM#2 (10 Mbps)  CCD Data
-------------------------
IDL>  file2=''   or  file2 = dialog_pickfile(filter='36*.*')  ; /Volumes/... for LOG file
IDL>  read_tm2, file2, /cd, /amegs, /bmegs

Enter time range (relative to T+0 in sec)  -  e.g.  0, 585

Quicklook:  Play the movies after images are extracted
----------
IDL>  movie_raw_megs, MEGS_FILE_NAME, 'A'

Note that MOVIE will not work if using IDL version 6.3 on Mac (IDL bug - get IDL v7.0).

Read Specific Image from Movie:
IDL>  aimage = 8  ; index into images in the file
IDL>  movie_raw_megs, MEGS_FILE_NAME, 'A', image=aimage
;  on return, "aimage" is the 8th image


TM#1 (3 Mbps)  Analog & Serial Data
-----------------------------------
IDL>  file1=''   or  file1 = dialog_pickfile(filter='36*.*')  ; /Volumes/... for DVD
IDL>  read_tm1, file1, /cd, /analog, /esp, /pmegs, /sps_xrs, /sps_owls, /cmd

Enter time range (relative to T+0 in sec)  -  e.g.  0, 585

Quicklook:  Plot the data after they are extracted
---------
IDL>  plot_analogs, ANALOG_FILE_NAME

IDL>  plot_esp, ESP_FILE_NAME

IDL>  plot_megsp, MEGSP_FILE_NAME

IDL>  plot_dualsps, channel, DUAL_SPS_FILE_NAME
					Channel names can be 'SPS', 'SPS_X', 'SPS_Y' and 'EUVx' where x=1-4 and X55
					There two DUAL SPS units on 36.353: XRS-X55 and OWLS-EUVOP



**** END OF FILE ****




************************************************************************************************
************************ EVERYTHING BELOW HERE IS FROM 36.336 2018 FLIGHT ***********************
************************************************************************************************

36.336 Launch Date = June 18, 2018 (DOY 169)  Time = 19:00:00 UT (68400.0 sec)
Apogee of 293.3 km at T+279.15 sec  (Radar Values)

Earth-Sun Distance = 1.01605 AU  (Earth Irradiance Correction to 1 AU = 0.9687)
F10.7 = 76.3  <F10.7> = 73.4   (MSIS: previous day F10.7=73.2, previous day Ap=5)
Ap = 20 (high activity)

IDL SETUP
---------
The alias rocket_idl, rocket_data need to be defined before running IDL.
Example definitions that can live in your .cshrc startup file.

alias rocket_idl 'setenv IDL_PATH "+/Users/twoods/Dropbox/minxss_dropbox/code/rocket_real_time:+/Users/twoods/Dropbox/minxss_dropbox/code/minxss_library:<IDL_DEFAULT>"'

setenv rocket_data /Users/twoods/Dropbox/minxss_dropbox/rocket_eve_docs/36.336/TM_Data

TM Flight Data Files are stored in $rocket_data/Flight/TM1 and $rocket_data/Flight/TM2

FLIGHT FILES FOR ROCKET INSTRUMENTS
------------
SDO EVE MEGS-A:	$rocket_data/Flight/TM2/flight_TM2_0_600_image_amegs.sav   (includes SAM)
SDO EVE MEGS-B:	$rocket_data/Flight/TM2/flight_TM2_0_600_image_bmegs.sav

SDO EVE ESP:	$rocket_data/Flight/TM1/flight_TM1_esp.sav

SDO EVE MEGS-P:	$rocket_data/Flight/TM1/flight_TM1_pmegs.sav

Rocket XRS:		$rocket_data/Flight/TM1/flight_TM1_goes_xrs.sav  (includes X123, SPS, PS-NIR)

CSOL Realtime:	$rocket_data/Flight/TM2/flight_TM2.log_0_600_csol.dat   (binary file)
CSOL 5-sec:		This appeared to fail (probably due to uplink issues during flight).
CSOL Full Images:  TBD file
					(SD-Card address = 200 to 247 for T+0 to T+600)
CSOL SPS:		$rocket_data/Flight/TM1/flight_TM1_sps_csol.txt

FPGA Commands:	$rocket_data/Flight/TM1/flight_TM1_cmd_fpga.txt

Analog Monitors: $rocket_data/Flight/TM1/flight_TM1_analogs.sav

FLIGHT RESULTS
--------------
There was a 30 Roll at T+350sec.  It did not work well so data after T+350s should be avoided.
To keep atmospheric absorption to less than 10%, the time period from T+170 to T+340 is best for analysis.

ESP Quad Diode results for flight near apogee:
    ESP X offset = +8.31 arc-mintues (wavelength shift axis for ESP) (SPARCS +Yaw)
    ESP Y offset = -4.78 arc-minutes  (SPARCS +Pitch)
    	NOTE: 3 active regions on the Sun during this flight (X-ray can skew center !)

SAM Center results:  (with ESP grating installed in SAM)
	SAM Sun Center X = 1573 pixel +/- 10 pixels    VERY DIFFICULT TO SEE SAM SUN CENTER
	               Y = 234 pixel  +/- 10 pixels
	SAM Perfect Center X = 1610, Y = 277
	SAM Angle (arc-min) = (pixel_offset * 15 microns / (32E4 microns))  * 180. * 60. / !PI
	SAM Offset X = -6.0 arc-min  +/- 2 arc-min
	           Y = -6.9 arc-min  +/- 2 arc-min

CSOL SPS Center results:
	Quad X = -1.2866    (post-vibe alignment: -1.293)   Shift of 0.38 arc-min
	Quad Y =  0.5563	(post-vibe alignment:  0.545)   Shift of 0.68 arc-min
	Quad Sum = 51773.

CSOL PicoSIM VIS are all saturated: all 6 channels went to 65535. then down to 31745.
	for rest of solar exposure.

X123 SPS Center results:
	Quad X =  0.421    (post-vibe alignment:  0.450)  Shift of -1.7 arc-min
	Quad Y = -0.232	   (post-vibe alignment: -0.190)  Shift of -2.5 arc-min
	Quad Sum = 55870.


Dark Measurements are at T+60 and T+490

Apogee at T+276 at about 280 km.

XRS/ESP Mechanism Closed at T+238 +/- 6 sec.

**************************   PROCESSING / PLOTTING CODE IN IDL  ********************************
There are two types of raw binary data files:
(1)  DVD (CD) files created as post-processing by WSMR  (complete data set)
(2)  DataView files created real-time (RT) (incomplete - not all data captured this way)

Instructions are different based on which file you are trying to process to extract out data.

Tom has processed these files (using instructions below) and extracted the data and saved into
*.dat binary files and also as *.sav IDL savesets for each instrument.
Others do not need to reprocess the Binary files.

NOTE:
IDL save set (*.sav) files are better to use than the binary data files due to new binary file
issues between different type computers in how they consider integers as either 32-bit or 64-bit integers.

Projects phase_development server: Rocket_Woods/Flights/36.290/Flight_Data/
and
EVESCI2 Directory:  /eve_analysis/testing/analysis/Rocket_flight_analysis/rocket_36290/

Subdirectories include:
	binary
	code
	pictures
	saveset
	video

As end user, see the following QUICKLOOK instructions for the data (*.dat) files.

**********   WARNING: while most code has been updated for 36.290 processing, there still might
be some code that works with the old rocket flights but not yet updated for 36.290 (2013 flight).

-----------------------------------------------------------
QUICKLOOK Instructions - How to plot / read the rocket data
-----------------------------------------------------------

MEGS CCD Images - movie & read single image
-------------------------------------------
IDL>  movie_raw_megs, MEGS_FILE_NAME, 'A' or 'B', info=info

Note that MOVIE will not work if using IDL version 6.3 on Mac (IDL bug - get IDL v7.0).

This procedure is written to be able to read single image at a time:
IDL>  aimage = 8  ; index into images in the file
IDL>  movie_raw_megs, MEGS_FILE_NAME, 'A', image=aimage
;  on return, "aimage" is the 8th image

MEGS_FILE_NAME = 36290_Flight_TM2_TTT_TTT_raw_Xmegs.sav

where        X = "a" or "b"
and
       TTT-TTT = -200_100   for dark and flat field data - 29 images
               = Ê100_460   for primary solar data (near apogee) - 34 images
               = Ê460-585   for more dark and flat field data at end of flight - 12 images


ESP Data - plot & read data
---------------------------
IDL>  plot_esp, '36290_Flight_TM1_0_585_esp.SAV', esp_data

MEGS-P Data - plot & read data
---------------------------
IDL>  plot_megsp, '36290_Flight_TM1_0_585_pmegs.SAV', pmegs_data

MEGS-P packet also has the CCD analog monitors.  For plotting those in addition to the
photometer data, use the /all option.

IDL>  plot_megsp, '36290_Flight_TM1_0_585_pmegs.SAV', pmegs_data, /all


XPS Data - plot & read data
---------------------------
IDL>  plot_xps, '36286-tm1-Flt_-200_585_xps.sav', xps_data


GOES XRS Data - plot & read data
---------------------------
IDL>  plot_flt_rxrs, '36290_Flight_TM1_0_585_goes_xrs.dat', data=data


Analog Monitors Data - plot & read data
----------------------------------------
IDL>  plot_analogs, '36290_Flight_TM1_0_585_analogs.SAV', analog_data



---------------------------------------------------
1.  Instructions for Extracting Data from DVD (CD) Binary Files (complete data set)
---------------------------------------------------

DVD TM File Names
-----------------
*.bin are TM files

36290_Flight_TM1.log  is TM#1 file (5 Mbps) for Analog monitors, ESP, XPS, MEGS-P, GOES-XRS

36290_Flight_TM2.log is TM#2 file (10 Mbps) for CCD data (MEGS-A & B)

TM#2 (10 Mbps)  CCD Data
-------------------------
IDL>  file2=''   or  file2 = dialog_pickfile(filter='36*.*')  ; /Volumes/... for DVD
IDL>  read_tm2, file2, /cd, /amegs, /bmegs, /csol

Enter time range (relative to T+0 in sec)  -  e.g.  0, 585

Quicklook:  Play the movies after images are extracted
----------
IDL>  movie_raw_megs, MEGS_FILE_NAME, 'A'

Note that MOVIE will not work if using IDL version 6.3 on Mac (IDL bug - get IDL v7.0).

Read Specific Image from Movie:
IDL>  aimage = 8  ; index into images in the file
IDL>  movie_raw_megs, MEGS_FILE_NAME, 'A', image=aimage
;  on return, "aimage" is the 8th image

Quicklook for CSOL:  Play the movies after images are extracted
-------------------
IDL>  movie_csol, CSOL_FILE_NAME


TM#1 (5 Mbps)  Analog & Serial Data
-----------------------------------
IDL>  file1=''   or  file1 = dialog_pickfile(filter='36*.*')  ; /Volumes/... for DVD
IDL>  read_tm1, file1, /cd, /analog, /esp, /pmegs, /xrs, /cmd, /sps_csol

Enter time range (relative to T+0 in sec)  -  e.g.  0, 585

Quicklook:  Plot the data after they are extracted
---------
IDL>  plot_analogs, ANALOG_FILE_NAME

IDL>  plot_esp, ESP_FILE_NAME

IDL>  plot_megsp, MEGSP_FILE_NAME

IDL>  plot_rxrs, GOES_FILE_NAME

IDL>  plot_picosim, channel, PICOSIM_SPS_FILE_NAME
					Channel names can be 'SPS', 'SPS_X', 'SPS_Y' and 'PSx' where x=1-6


---------------------------------------------------
2.  Instructions for  Extracting Data from DataView (RT) Files (incomplete data set)
---------------------------------------------------

DataView TM File Names
-----------------
Raw_Data_TM1_YY_MM_DD_HH-MM* are files for TM#1 (analog monitors, ESP, XPS, MEGS-P)

Raw_Data_TM2_YY_MM_DD_HH-MM* are files for TM#2 (CCD data for MEGS-A & B)


TM#2 (10 Mbps)  CCD Data
-------------------------
IDL>  file2=''   or  file2 = dialog_pickfile(filter='Raw_Data_TM2*')
IDL>  read_tm2, file2, /amegs, /bmegs, /csol

Enter time range (relative to T+0 in sec)  -  e.g.  0, 450

Quicklook:  Play the movies after images are extracted
----------
IDL>  movie_megs, MEGS_FILE_NAME, 'A', info=ainfo

Note that MOVIE will not work if using IDL version 6.3 on Mac (IDL bug - get IDL v7.0).

Read Specific Image from Movie:
IDL>  aimage = 8  ; index into images in the file
IDL>  movie_megs, MEGS_FILE_NAME, 'A', image=aimage
;  on return, "aimage" is the 8th image


TM#1 (5 Mbps)  Analog & Serial Data
-----------------------------------
IDL>  file1=''   or  file1 = dialog_pickfile(filter='Raw_Data_TM1*')
IDL>  read_tm1, file1, /analog, /esp, /pmegs, /xrs, /cmd, /sps_csol

Enter time range (relative to T+0 in sec)  -  e.g.  0, 450

Quicklook:  Plot the data after they are extracted
---------
IDL>  plot_analogs, ANALOG_FILE_NAME

IDL>  plot_esp, ESP_FILE_NAME

IDL>  plot_megsp, MEGSP_FILE_NAME

IDL>  plot_rxrs, GOES_FILE_NAME

IDL>  plot_picosim, channel, PICOSIM_SPS_FILE_NAME
					Channel names can be 'SPS', 'SPS_X', 'SPS_Y' and 'PSx' where x=1-6

**** END OF FILE ****


