# test ModelSetup/Medium/mesh.jl
include("typedefine.jl")
include("mesh.jl")

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

#=== Case 2 layered medium ===#
vp = 5000:100:6000
vs = 3500:100:4500
rho = 2.8
depth = 100:100:1100
Horizon = 1000
pkf = 10
T = 2
ext = 10
iflag = 2
dx = nothing
dz = nothing
dt = nothing
Tn = nothing

Model = model(vp,vs,rho,depth,Horizon,pkf,T,ext,iflag,dx,dz,dt,Tn)

#=== Case 3 homogeneous medium ===#
vp = 5000
vs = 3500
rho = 2.8
depth = 1000
Horizon = 1000
pkf = 10
T = 2
ext = 10
iflag = 2
dx = nothing
dz = nothing
dt = nothing
Tn = nothing

Model = model(vp,vs,rho,depth,Horizon,pkf,T,ext,iflag,dx,dz,dt,Tn)
