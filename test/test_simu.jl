using seisimu

vp = [5000.0]
vs = [3500.0]
rho = 2.8
depth = 3000
Horizon = 3000
pkf = 30
T = 7
ext = 10
iflag = 2
dx = nothing
dz = nothing
dt = nothing
nT = nothing

model = initmodel(vp,rho,depth,Horizon,pkf,T)

sn = 1
loc = [1000 1000]
ot = [0.001]
mt = [1 1 0]
wavelet = Ricker(pkf,model.medium.dt)

sou = initsource(sn,loc,ot,mt,wavelet,model)

@time runsimu(model,sou)
