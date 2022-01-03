# TO RUN:
# include("main.jl")
#
# AUTHOR
# Dr. Umberto Saetti, Assistant Professor, Department of Aerospace Engineering,
# Auburn university
# saetti@auburn.edu
#
# LAST UPDATED
# 1/3/2022
#
# DESCRIPTION
# Nonlinear simulation model of the rigid-body dynamics of a conventional
# helicopter configuration (i.e., main rotor + tail rotor).
# The model uses static aerodynamic models for fuselage, main rotor, tail
# rotor, and empennage. The trim, linearization routines, and overall
# architecture of the simulation follows the teachings of Dr. Joe Horn at Penn
# State. A more detailed description of the model, including references can be
# found in Ref. 1. If using the code, please cite Ref. 1. 
#
# REFERENCES
# 1) Saetti, U., and Horn, J. F., "Flight Simulation and Control using the Julia
#    Language", AIAA Scitech Forum, San Diego, CA, Jan 3-7, 2022.
#    DOI: https://arc.aiaa.org/doi/10.2514/6.2022-2354.
#
# STATES             | NINDICES    | DESCRIPTION
# ______________________________________________________________________________
# u v w              |  1  2  3    | body velocities [ft/s]
# p q r              |  4  5  6    | angula rates [rad]
# phi theta psi      |  7  8  9    | Euler angles [rad]
# x  y  z            | 10 11 12    | position [ft]
#
# CONTROL INPUTS     | INDICED     | DESCRIPTION
# ______________________________________________________________________________
# delta_lat          | 1           | lateral stick [%]
# delta_lon          | 2           | longitudinal stick [%]
# delta_col          | 3           | collective stick [%]
# delta_ped          | 4           | pedals [%]
#
# ------------------------------------------------------------------------------

# include packages
#using DifferentialEquations
using LinearAlgebra
using SparseArrays
using ControlSystems
using MAT
using Interpolations
using ApproxFun
using Printf
#ENV["MPLBACKEND"]="tkagg"
using Plots
using JLD
#using PyPlot
#using GR

@printf("\n                 J-SimpleHel           \n")
@printf("\n          Author: Dr. Umberto Saetti\n")
@printf("\n-----------------------------------------------\n")


# include functions
include("atan2.jl")
include("linearize.jl")
include("trimmer.jl")
include("SimpleHel.jl")
include("control_mixing.jl")
include("SimpleHel_DI.jl")
include("rk4.jl")
include("finp.jl")
include("simulate.jl")
include("simulate_DI.jl")
include("wrapper.jl")
# load H-60 constants
include("H60_constants.jl")

# ------------------------ TRIM FLIGHT DYNAMICS MODEL --------------------------

# trim speeds [ft/s]
VXTRIM=80*1.688
VYTRIM=0.0
VZTRIM=0.0
# trim turn rate [rad/s]
PSIDTRIM=0.0
include("trimInit.jl")
# trim aircraft
x0, u0 = trimmer(SimpleHel!,x0,u0,xdot_targ)
# linearize aircraft
A, B, C, D = linearize(SimpleHel!,x0,u0)
# partition in longitudinal and lateral flight dynamics
Alon=A[[1,3,5,8],[1,3,5,8]]
Alat=A[[2,4,6,7],[2,4,6,7]]

# ----------------------------- SPECTRAL ANALYSIS ------------------------------

# eigenvalues
eigs=eigvals(A)
eigsLon=eigvals(Alon)
eigsLat=eigvals(Alat)
# store eigenvalues
#save("AB_J-SimpleHel.jld","A",A,"B",B)
#save("eigs_J-SimpleHel.jld","eigs",eigs)

# plot eigenvalues
gr(size=(800,600))
plt_eigs=plot(real(eigs),imag(eigs),seriestype=:scatter,
      label="Full-Order", legendfontsize=12, legend=:topleft,
      marker = (:circle, :royalblue, 6), markerstrokecolor = :royalblue)
plot!(real(eigsLat),imag(eigsLat),seriestype=:scatter,
      label="Reduced-Order (Lat)", legendfontsize=12,
      marker = (:utriangle, :brown3, 6), markerstrokecolor = :brown3,
      reuse = false, color=:purple)
plot!(real(eigsLon),imag(eigsLon),seriestype=:scatter,
      label="Reduced-Order (Lon)", legendfontsize=12,
      marker = (:star5, :forestgreen, 6), markerstrokecolor = :forestgreen)
display(plt_eigs)
xaxis!("Real",xguidefontsize=14)
yaxis!("Imag",yguidefontsize=14)

# ---------------------------- FREQUENCY RESPONSES -----------------------------

# linear system
sys=ss(A,B,Matrix{Float64}(I, NSTATES, NSTATES),zeros(NSTATES,NCTRLS))
# transfer function (p/dlat)
pdlat=tf(sys[4,1])
# freq response
mag, phase, w = bode(pdlat, 10.0.^(range(-1,stop=2,length=1000)))

# plot pdlat
gr(size=(1000,1000))
plt_mag=plot(vec(w), 20*log10.(vec(mag)/perc2in[1,1]), label="Full order",
             reuse = false, line=(4, :black, :solid))
xaxis!(:log10)
yaxis!("Magnitude [dB]")
if VXTRIM==80*1.688
      title!("p/dlat, KTAS = 80 kts")
elseif VXTRIM<=1*1.688
      title!("p/dlat, KTAS = 0 kts")
end
xlims!((0.1,100))
plt_phase=plot(vec(w), wrapper(vec(phase)).-360, label="", line=(4, :black, :solid))
yaxis!("Phase [deg]")
xaxis!(:log10)
xlims!((0.1,100))
plt_pdlat=plot(plt_mag, plt_phase, layout=(2,1), margin=10Plots.mm)
display(plt_pdlat)
#if VXTRIM==80*1.688
#      savefig(".\\Plots\\pdlat_80kt.png")
#elseif VXTRIM<=1*1.688
#      savefig(".\\Plots\\pdlat_Hover.png")
#end

# transfer function (q/dlon)
pdlat=tf(sys[5,2])
# freq response
mag, phase, w = bode(pdlat, 10.0.^(range(-1,stop=2,length=1000)))

# plot qdlon
gr(size=(1000,1000))
plt_mag=plot(vec(w), 20*log10.(vec(mag)/perc2in[2,2]), label="Full order",
             reuse = false, line=(4, :black, :solid))
xaxis!(:log10)
yaxis!("Magnitude [dB]")
if VXTRIM==80*1.688
      title!("q/dlon, KTAS = 80 kts")
elseif VXTRIM<=1*1.688
      title!("q/dlon, KTAS = 0 kts")
end
xlims!((0.1,100))
plt_phase=plot(vec(w), wrapper(vec(phase)).-360, label="", line=(4, :black, :solid))
yaxis!("Phase [deg]")
xaxis!(:log10)
xlims!((0.1,100))
plt_qdlon=plot(plt_mag, plt_phase, layout=(2,1), margin=10Plots.mm)
display(plt_qdlon)
#if VXTRIM==80*1.688
#      savefig(".\\Plots\\qdlon_80kt.png")
#elseif VXTRIM<=1*1.688
#      savefig(".\\Plots\\qdlon_Hover.png")
#end

# ------------------------------ TIME SIMULATION -------------------------------

# time step [s]
dt=0.01
# length of simulation [s]
Tsim=10
# run open-loop simulation
state_OL, time_OL = simulate(SimpleHel!,finp,Tsim,dt,x0,u0)
# number of dynamic inverse states
NDISTATES=9
# initial state vector of closed-loop simulation
x0_DI=[x0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; u0[1]; u0[2]; u0[4]]
# run closed-loop simulation
state_CL, time_CL = simulate_DI(SimpleHel_DI!,finp,Tsim,dt,x0,u0)
# store open- and closed-loop simulations
#save("resp_SimpleHel.jld","time_OL",time_OL,"time_CL",time_CL,
#    "state_OL",state_OL,"state_CL",state_CL)

# plot attitude
gr(size=(800,600))
i=7
p1=plot(time_OL,state_OL[i,:]*180/pi,label="Open Loop",legend=:topright, legendfontsize=12,line=(3, :royalblue, :solid))
plot!(time_OL,state_CL[i,:]*180/pi,label="Closed Loop",legend=:topright, legendfontsize=12,line=(3, :brown3, :dash))
yaxis!("φ [deg]",yguidefontsize=14)
xlims!((0,10))
i=8
p2=plot(time_OL,state_OL[i,:]*180/pi,label="",line=(3, :royalblue, :solid))
plot!(time_OL,state_CL[i,:]*180/pi,label="",line=(3, :brown3, :dash))
yaxis!("θ [deg]",yguidefontsize=14)
xlims!((0,10))
i=9
p3=plot(time_OL,state_OL[i,:]*180/pi,label="",line=(3, :royalblue, :solid))
plot!(time_OL,state_CL[i,:]*180/pi,label="",line=(3, :brown3, :dash))
xaxis!("Time [s]",xguidefontsize=14)
yaxis!("ψ [deg]",yguidefontsize=14)
xlims!((0,10))
p4=plot(p1, p2, p3, layout=(3,1), leftmargin=3Plots.mm)
display(p4)
#savefig(".\\Plots\\att_J-SimpleHel.svg")

# plot angular rates
gr(size=(800,600))
i=4
p5=plot(time_OL,state_OL[i,:],label="Open Loop",legend=:bottomleft, legendfontsize=12,line=(3, :royalblue, :solid))
plot!(time_OL,state_CL[i,:],label="Closed Loop",legend=:bottomleft, legendfontsize=12,line=(3, :brown3, :dash))
yaxis!("p [rad/s]",yguidefontsize=14)
xlims!((0,10))
i=5
p6=plot(time_OL,state_OL[i,:],label="",line=(3, :royalblue, :solid))
plot!(time_OL,state_CL[i,:],label="",line=(3, :brown3, :dash))
yaxis!("q [rad/s]",yguidefontsize=14)
xlims!((0,10))
i=6
p7=plot(time_OL,state_OL[i,:],label="",line=(3, :royalblue, :solid))
plot!(time_OL,state_CL[i,:],label="",line=(3, :brown3, :dash))
xaxis!("Time [s]",xguidefontsize=14)
yaxis!("r [rad/s]",yguidefontsize=14)
xlims!((0,10))
p8=plot(p5, p6, p7, layout=(3,1), leftmargin=3Plots.mm)
display(p8)
#savefig(".\\Plots\\ang_J-SimpleHel.svg")
