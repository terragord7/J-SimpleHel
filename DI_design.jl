# To run:
# include("DI_design.jl")

# DESCRIPTION
# Trims and linearizes the rotorcraft model at incremental flight speeds, and
# generates the gains of the dynamic inversion controller. Always run after
# modifying the aircraft parameters in the AIRCRAFTNAME_constants.jl file.

# include packages
#using DifferentialEquations
using LinearAlgebra
using SparseArrays
using ControlSystems
using MAT
using Interpolations
using ApproxFun
using Printf
using Plots
using JLD
using LaTeXStrings

# include functions
include("atan2.jl")
include("linearize.jl")
include("trimmer.jl")
include("SimpleHel.jl")
#include("rk4.jl")
#include("finp.jl")
# load H-60 constants
include("H60_constants.jl")
# vector of trim speeds [kts]
VXTRIM_vec=[0.01,20.0,40.0,60.0,80.0,100.0,120.0,140.0,160.0]
# initialize DI matrices
CATAB = zeros(3,3,length(VXTRIM_vec))
CBinvTAB = zeros(3,3,length(VXTRIM_vec))
x0_mat = zeros(NSTATES,length(VXTRIM_vec))
u0_mat = zeros(NCTRLS,length(VXTRIM_vec))
# trim at different speeds
for iv=1:length(VXTRIM_vec)
    # load H60 constants
    include("H60_constants.jl")
    # forward, lateral, and velocities in trim [ft/s]
    VXTRIM=VXTRIM_vec[iv]*1.688
    VYTRIM=0.0
    VZTRIM=0.0
    # turn rate in trim [rad/s]
    PSIDTRIM=0.0
    # initial heading [rad]
    PSITRIM = 0.0*pi/180
    # initial position [ft]
    ALTTRIM = 0.0
    XNTRIM = 0.0
    YETRIM = 0.0
    # initialize state vectors and initial guess
    x0=zeros(NSTATES)
    x0[1]=VXTRIM
    x0[2]=VYTRIM
    x0[3]=VZTRIM
    # option to trim in zero sideselip or zero bank angle
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
    # time [s]
    t=0.0
    # trim aircraft
    x0, u0 = trimmer(SimpleHel!,x0,u0,xdot_targ)
    # linearized aircraft dynamics
    A, B, C, D = linearize(SimpleHel!,x0,u0)
    # store trim state and controls
    x0_mat[:,iv] = x0
    u0_mat[:,iv] = u0
    # dynamic inverse matrices
    C_DI=[1 0 0
          0 1 0
          0 0 1]
    CATAB[:,:,iv]=C_DI*A[4:6,4:6]
    CBinvTAB[:,:,iv]=inv(C_DI*B[4:6,[1,2,4]])
end

# command filters time constants
tau_p=1/3.5
tau_q=1/4.5
tau_r=0.5
# washout filter for controls
tau_u=5.0;
# disturbance rejection frequency and damping and integrator poles
wnroll_d=3.5
dmproll_d=1.0
wnpitch_d=4.5
dmppitch_d=1.0
wnyaw_d=2.;
dmpyaw_d=1.;
# feedback gains
kp_roll=2*wnroll_d*dmproll_d
ki_roll=wnroll_d^2
kp_pitch=2*wnpitch_d*dmppitch_d
ki_pitch=wnpitch_d^2
kp_yaw=2*dmpyaw_d*wnyaw_d
ki_yaw=wnyaw_d^2;
# pilot gains
k_lat=pi/100
k_lon=pi/100
k_ped=0.644
# save DI matrices and gains into a file
save("DI.jld","CATAB",CATAB,"CBinvTAB",CBinvTAB,"kp_roll",kp_roll,"ki_roll",
    ki_roll,"kp_pitch",kp_pitch,"ki_pitch",ki_pitch,"kp_yaw",kp_yaw,"ki_yaw",
    ki_yaw,"k_lat",k_lat,"k_lon",k_lon,"k_ped",k_ped,"tau_p",tau_p,
    "tau_q",tau_q,"tau_r",tau_r,"tau_u",tau_u)

# store trim state and control input vectors
save("trim_J-SimpleHel.jld","x0_mat",x0_mat,"u0_mat",u0_mat,"VXTRIM_vec",VXTRIM_vec)

# plot trim attitudes
#gr(size=(1000,1000))
gr(size=(800,600))
p1=plot(VXTRIM_vec,x0_mat[7,:]*180/pi,label="J-GenHel",legend=:bottomright, legendfontsize=12,line=(3, :blue, :solid), marker = (:circ, :blue, 8), markerstrokecolor = :blue)
xaxis!((0, 160), 0:20:160)
yaxis!("φ [deg]",yguidefontsize=14)
p2=plot(VXTRIM_vec,x0_mat[8,:]*180/pi,label="",line=(3, :blue, :solid), marker = (:circ, :blue, 8), markerstrokecolor = :blue)
xaxis!((0, 160), 0:20:160)
yaxis!("θ [deg]",yguidefontsize=14)
p3=plot(VXTRIM_vec,x0_mat[9,:]*180/pi,label="",line=(3, :blue, :solid), marker = (:circ, :blue, 8), markerstrokecolor = :blue)
xaxis!((0, 160), 0:20:160)
xaxis!("Total speed, V [kts]",xguidefontsize=14)
yaxis!("ψ [deg]",yguidefontsize=14)
p4=plot(p1, p2, p3, layout=(3,1))
display(p4)
#savefig(".\\Plots\\trim.svg")

# plot trim controls
gr(size=(800,600))
p5=plot(VXTRIM_vec,u0_mat[1,:],label="J-SimpleHel",legend=:bottomright, legendfontsize=12,line=(3, :blue, :solid), marker = (:circ, :blue, 8), markerstrokecolor = :blue)
xaxis!((0, 160), 0:20:160)
yaxis!(L"\delta_\mathrm{lat} [\%]",yguidefontsize=14)
p6=plot(VXTRIM_vec,u0_mat[2,:],label="",line=(3, :blue, :solid), marker = (:circ, :blue, 8), markerstrokecolor = :blue)
xaxis!((0, 160), 0:20:160)
yaxis!(L"\delta_\mathrm{lon} [\%]",yguidefontsize=14)
p7=plot(VXTRIM_vec,u0_mat[3,:],label="",line=(3, :blue, :solid), marker = (:circ, :blue, 8), markerstrokecolor = :blue)
xaxis!((0, 160), 0:20:160)
yaxis!(L"\delta_\mathrm{col} [\%]",yguidefontsize=14)
p8=plot(VXTRIM_vec,u0_mat[4,:],label="",line=(3, :blue, :solid), marker = (:circ, :blue, 8), markerstrokecolor = :blue)
xaxis!((0, 160), 0:20:160)
xaxis!("Total speed, V [kts]",xguidefontsize=14)
yaxis!(L"\delta_\mathrm{ped} [\%]",yguidefontsize=14)
p9=plot(p5, p6, p7, p8, layout=(4,1))
display(p9)
#savefig(".\\Plots\\trim_ctrls.svg")
