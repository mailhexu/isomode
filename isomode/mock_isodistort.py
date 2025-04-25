"""Mock ISODISTORT service for testing"""

SAMPLE_MODE_DETAIL = """Parent structure (123 P4/mmm)
a=3.92640, b=3.92640, c=15.18110, alpha=90.00000, beta=90.00000, gamma=90.00000

Displacive mode definitions
P4/mmm[0,0,0]GM1+(a)[Pb1:h:dsp]A1(a) normfactor = 0.01647
P4/mmm[0,0,0]GM5+(a,b)[Pb1:h:dsp]E(a) normfactor = 0.12734
P4/mmm[0,0,0]GM5-(a,b)[K1:d:dsp]Eu(a) normfactor = 0.09005

Displacive mode amplitudes
mode                               As        Ap       dmax
[0,0,0]GM1+  all                  0.00052  0.00037
[0,0,0]GM5+  all                  0.00000  0.00000
[0,0,0]GM5-  all                  0.00000  0.00000
"""

def get_mode_details():
    """Return sample mode detail text"""
    return SAMPLE_MODE_DETAIL
