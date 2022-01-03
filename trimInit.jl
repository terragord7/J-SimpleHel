# DESCRIPTION
# Trim initialization.
#
# ------------------------------------------------------------------------------

# initial guess for state in trim solution
x0=zeros(12,1)
x0[1]=VXTRIM
x0[2]=VYTRIM
x0[3]=VZTRIM
if (VXTRIM>=60*1.688)
    if PSIDTRIM==0
        # Trim heading for zero bank angle above 60 kts
        x0[7]=0.0;
        TRIMVARS[7:8]=[8; 9]
    else
        x0[9]=atan2(VYTRIM,VXTRIM)
        TRIMVARS[7:8]=[7; 8]
        # initial guess for phi (using turn coord)
        x0[7]=atan(PSIDTRIM*sqrt(VXTRIM^2+VYTRIM^2+VZTRIM^2)/G)
        # initial guess for r
        x0[6]=G/sqrt(VXTRIM^2+VYTRIM^2+VZTRIM^2)*sin(x0[7])
        # initial guess for q
        x0[5]=x0[6]*tan(x0[7])
    end
else
    # otherwise set heading to flight path
    x0[9]=atan2(VYTRIM,VXTRIM)
    TRIMVARS[7:8]=[7; 8]
end
# initial guess for controls [deg]
u0=[50.0; 50.0; 50.0; 50.0]
# target values for state derivatives
xdot_targ=zeros(12)
xdot_targ[9]=PSIDTRIM
xdot_targ[10]=VXTRIM
xdot_targ[11]=VYTRIM
xdot_targ[12]=VZTRIM
