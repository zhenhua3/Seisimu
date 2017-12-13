using seisimu

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

model = initmodel(vp,vs,rho,depth,Horizon,pkf,T)
