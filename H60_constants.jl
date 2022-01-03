# DESCRIPTION
# This file contains the helicopter parameters.. Always re-run DI_design.jl
# after modifying the aircraft parameters.

# ------------------------- MISCELLANEOUS CONSTANTS ----------------------------

# air density [slugs-ft³]
RHO=0.002378
# acceleration of gravity [ft/s²]
g=32.17

# ----------------------------- MASS PROPERTIES --------------------------------

# gross weight [lb]
W=16000.0
# inertias [sl-ft²]
Ix=4659.0*W/16000.0
Iy=38512.0*W/16000.0
Iz=36796.0*W/16000.0
Ixz=1882.0*W/16000.0
# CG location (station, water line, butt line) [in]
FSCG=358.0
BLCG=0.0
WLCG=248.2

# --------------------------------- FUSELAGE -----------------------------------

# fuselage drag [ft²]
CDSF=35.04
# location of fuselage aerodynamic center (station, water line, butt line) [in]
FSfus=345.5
BLfus=0.0
WLfus=234.0
# number of blades
NB=4;
# single blade weight [lb]
WBLADE=256.9
# aircraft gross weight (rotor blades not included) [lb]
WEIGHTNR = W-NB*WBLADE;
# aircraft mass (rotor blades not included) [slug]
mass=WEIGHTNR/g;
# main rotor location (station, water line, butt line) [in]
FSMR=341.25
BLMR=0.0
WLMR=315.0
# fuselage CG location (without rotor) (station, water line, butt line) [in]
FSCGB=((W*FSCG)-(NB*WBLADE*FSMR))/WEIGHTNR
WLCGB=((W*WLCG)-(NB*WBLADE*WLMR))/WEIGHTNR
BLCGB=((W*BLCG)-(NB*WBLADE*BLMR))/WEIGHTNR
# moment arms from fuselage CG to fuselage aero center [ft]
FWT = (FSCG-FSfus)/12
WWT = (WLCG-WLfus)/12
BWT = (BLCG-BLfus)/12
# location of fuselage aero center relative to CG [ft]
rFS=[FWT; BWT; WWT];

# -------------------------------- MAIN ROTOR ----------------------------------

# rotor Radius [ft]
RADIUS=26.83
# rotor speed [rad/s]
OMEGA=27.0
# solidity
SOL=0.0821
# blade lift slope [1/rad]
A0=5.73
# blade twist [rad]
TWIST=-13.0*pi/180.0
# lock number = rho*a0*c*R^4/Ibeta
GAMMA=8.0755
# moment arms from fuselage CG to main rotor hub [ft]
Xh=(FSCGB-FSMR)/12.0
Yh=(BLCGB-BLMR)/12.0
Zh=(WLCGB-WLMR)/12.0
# main rotor location
rMR=[Xh; Yh; Zh]
# rotor blade profile CD parameters
DELTA0=0.01
DELTA2=250.0

# -------------------------------- TAIL ROTOR ----------------------------------

# tail rotor speed [rad/s]
OMEGTR=124.62
# tail rotor radius [ft]
RTR=5.5
# tail rotor blade lift slope [1/rad]
A0TR=5.73
# tail rotor solidity
SOLTR=0.1875
# tail rotor twist [rad]
TWISTTR=-17.2*pi/180.0
# location (station, water line, butt line) [in]
FSTR=732.0
BLTR=-14.0
WLTR=324.7
# moment arms from fuselage CG to tail rotor hub
XTR=-(FSTR-FSCGB)/12.0
YTR=(BLTR-BLCGB)/12.0
ZTR=-(WLTR-WLCGB)/12.0
# location of tail rotor relative to CG [ft]
rTR=[XTR; YTR; ZTR]

# -------------------------- HORIZONTAL STABILIZER -----------------------------

# area [ft²]
SHT=45.0
# lift slope [1/rad]
AHT=2.3
# horizontal stabilizer location (station, water line, butt line) [in]
FSHT=700.0
BLHT=0.0
WLHT=244.0
# moment arms from fuselage CG to horizontal stabilizer aero center
XHT=-(FSHT-FSCG)/12.0
YHT=(BLHT-BLCG)/12.0
ZHT=-(WLHT-WLCG)/12.0
# location of tail rotor relative to CG [ft]
rHT=[XHT; YHT; ZHT]

# --------------------------- VERTICAL STABILIZER ------------------------------

# area [ft²]
SVT=32.3
# lift slope [1/rad]
AVT=2.3
# horizontal stabilizer location (station, water line, butt line) [in]
FSVT=695.0
BLVT=0.0
WLVT=273.0
# moment arms from fuselage CG to vertical stabilizer aero center
XVT=-(FSVT-FSCG)/12.0
YVT=(BLVT-BLCG)/12.0
ZVT=-(WLVT-WLCG)/12.0
# location of tail rotor relative to CG [ft]
rVT=[XVT; YVT; ZVT]

# ------------------------------ CONTROL MIXING --------------------------------

# conversion from stick percentages to stick inches
perc2in = [0.1 0   0   0
           0   0.1 0   0
           0   0   0.1 0
           0   0   0   5.38/100]
# mixer gain matrix
LNKGAIN = [0.2062 0      0      0
           0      0.2172 0      0
           0      0      0.2025 0
           0      0      0      0.3780]
# convert stick inputs to servo commands
MIXGAIN = [0.6800  0      0.9930   0
           0       1.2020 0.9050  -0.4288
           0      -1.2020 1.2990   0.4288
           0       0      -0.8554  1.6043]
# gain for conversion servo to blade deflections
SWASHGAIN = [11.3413 -5.6706 -5.6706
              0      -5.6706  5.6706
              0       3.5034  3.5034]
TRGAIN = -9.142
# bias for conversion servo to blade deflections
SWASHBIAS = [-7.9738, 9.7280, 10.1111]
TRBIAS = 21.98

# --------------------- LINEARIZATION AND TRIM PARAMETERS ----------------------

# number of states
NSTATES=12
# number of control inputs
NCTRLS=4
# number of output variables
NOUT=6
# state vector perturbations for linearization
DELXLIN=[0.1; 0.1; 0.1; pi/180*0.1*ones(6); 0.1; 0.1; 0.1]
# control vector perturbations for linearization
DELCLIN=[0.1; 0.1; 0.1; 0.1]
# trim targets
TRIMTARG=collect(1:12)
# trim variables
TRIMVARS=[collect(1:8); collect(13:16)]
# tolerance on trim error
TOL=0.001*[0.001; 0.001; 0.001; pi/180*0.001*ones(6); 0.001; 0.001; 0.001;
           0.001; 0.001; 0.001]
