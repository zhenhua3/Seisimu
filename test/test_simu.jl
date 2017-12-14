using seisimu

vp = [5000.0]
vs = [3500.0]
rho = 2.8
depth = 10000
Horizon = 10000
pkf = 30
T = 7
ext = 10
iflag = 2
dx = nothing
dz = nothing
dt = nothing
nT = nothing

model = initmodel(vp,vs,rho,depth,Horizon,pkf,T)

sn = 1
loc = [5000 5000]
ot = [0.001]
mt = [1 5 1]
wavelet = Ricker(pkf,model.medium.dt)

sou = initsource(sn,loc,ot,mt,wavelet,model)

@time runsimu(model,sou; showevery = 20000)
