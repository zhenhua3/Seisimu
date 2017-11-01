using SeismicBox

vp = [3500.0,3200,4000,4500]
vs = nothing
rho = 2.8
depth = [300,200,100,300]
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


sn = 2
loc = [500 500; 800 300]
ot = [0.8, 1.5]
tp = ["expl", "DC"]
wavelet = Ricker(pkf,model.medium.dt)
plt.show(plt.plot(wavelet))

waveform = [wavelet wavelet]

sou = initsource(sn,loc,ot,tp,waveform,model)

plt.show(plt.plot(sou[1].waveform))
