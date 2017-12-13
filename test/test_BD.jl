# test ModelSetup/Medium/initwf.jl
include("typedefine.jl")
include("mesh.jl")
include("initwf.jl")
include("BD.jl")

# 2d elastic #
vp = [5000.0,5300,5800,6000]
vs = [3500.0,3200,4000,4500]
rho = 2.8
depth = [700,200,400,600]
HX = 1000
pkf = 10
T = 2
ext = 10
iflag = 2
dx = nothing
dz = nothing
dt = nothing
nT = nothing

Model = model(vp,vs,rho,depth,HX,pkf,T,ext,iflag,dx,dz,dt,nT)

elwf = initnelwf(Model.nDZ,Model.nHX,ext,iflag)

pml = BD(elwf,Model)

# 2d acoustic #
vp = [5000.0,5300,5800,6000]
rho = 2.8
depth = [700,200,400,600]
HX = 1000
pkf = 10
T = 2
ext = 10
iflag = 2
dx = nothing
dz = nothing
dt = nothing
nT = nothing

Model = model(vp,rho,depth,HX,pkf,T,ext,iflag,dx,dz,dt,nT)

acwf = initnacwf(Model.nDZ,Model.nHX,ext,iflag)

pml = BD(acwf,Model)

# 3d elastic #
vp = [5000.0,5300,5800,6000]
vs = [3500.0,3200,4000,4500]
rho = 2.8
depth = [700,200,400,600]
HX = 1000
HY = 1000
pkf = 10
T = 2
ext = 10
iflag = 2
dx = nothing
dy = nothing
dz = nothing
dt = nothing
nT = nothing

Model = model(vp,vs,rho,depth,HX,HY,pkf,T,ext,iflag,dx,dy,dz,dt,nT)

elwf = initnelwf(Model.nDZ,Model.nHX,Model.nHY,ext,iflag)

pml = BD(elwf,Model)

# 3d acoustic #
vp = [5000.0,5300,5800,6000]
rho = 2.8
depth = [700,200,400,600]
HX = 1000
HY = 1000
pkf = 10
T = 2
ext = 10
iflag = 2
dx = nothing
dy = nothing
dz = nothing
dt = nothing
nT = nothing

Model = model(vp,rho,depth,HX,HY,pkf,T,ext,iflag,dx,dy,dz,dt,nT)

elwf = initnacwf(Model.nDZ,Model.nHX,Model.nHY,ext,iflag)

pml = BD(elwf,Model)
