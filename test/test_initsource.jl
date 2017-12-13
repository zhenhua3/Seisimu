include("CreateMesh.jl")

#=== Case 1 layered medium ===#
vp = [5000.0,5300,5800,6000]
vs = [3500.0,3200,4000,4500]
rho = 2.8
depth = [700,200,400,600]
Horizon = 1000
pkf = 10
T = 2
ext = 10
iflag = 2
dx = nothing
dz = nothing
dt = nothing
nT = nothing

Model = model(vp,vs,rho,depth,Horizon,pkf,T,ext,iflag,dx,dz,dt,nT)

nwf = initnelwf(Model.nDZ,Model.nHX,Model.ext,Model.iflag)

sn = 2
loc = [20 20;30 30]
ot = [0.1, 0.2]
mt = zeros(2,3)
waveform = Ricker(30,Model.dt)
waveform = [waveform waveform]

sou = initsource(sn,loc,ot,mt,waveform,Model,nwf)
