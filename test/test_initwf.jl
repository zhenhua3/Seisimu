# test ModelSetup/Medium/initwf.jl
include("typedefine.jl")
include("initwf.jl")


nDZ = 100
nHX = 100
ext = 10
iflag = 2

elwf = initelwf(nDZ,nHX,ext,iflag)

acwf = initacwf(nDZ,nHX,ext,iflag)

nHY = 100

elwf3d = initelwf(nDZ,nHX,nHY,ext,iflag)

acwf3d = initacwf(nDZ,nHX,nHY,ext,iflag)
