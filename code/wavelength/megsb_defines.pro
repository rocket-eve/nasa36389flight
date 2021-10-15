;
; MEGSB_DEFINES.PRO
;

MEGSB_GRATING_D = 1.e6/600. ; grating "d" in nm (or nm/line)
MEGSB_ANODE_WIDTH = .001    ; anode width in cm (pixel width)
MEGSB_FOCAL_LEN = 25.0      ; grating focal length in cm

machine_f=machar()
machine_d=machar(/double)


OK_STATUS = 0L
BAD_PARAMS = -1L
BAD_FILE = -2
BAD_FILE_DATA = -3
FILE_ALREADY_OPENED = -4
BAD_FOV_ANGLE = -5
BAD_CAL_FIT = -6

BAD_DATA_DIVIDE_BY_ZERO = -7

NUM_MEGSB_PIXELS = 2048L
NUM_MEGSB_COLS = 2048L
