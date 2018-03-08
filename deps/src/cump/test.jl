using seisimu
using PyCall
@pyimport matplotlib.pyplot as plt
include("extrapcump.jl")
pvel = 5000;
rho = 2.8;
DZ = 4000;
HX = 4000;
HY = 100;
pkf = 10;
T = 10;
ext = 20;
iflag = 2;

# === cuda ===#
model1 = initmodel(pvel,nothing,rho,DZ,HX,HY,pkf,T);

CPU_mem = 14
GPU_mem = 2.8
num_stream = 2
threadim = [8,4,8]

sn = 1
loc = [2000 2000 50]
ot = [0.05]
mt = [1 1 1 0 0 0]
wavelet = Ricker(pkf,model1.medium.dt)
waveform = wavelet
sou = initsource(sn, loc, ot, mt, waveform, model1)

@time extrap!(model1, threadim, CPU_mem, GPU_mem; num_stream = num_stream)

#=== openmp ===#
model2 = initmodel(pvel,nothing,rho,DZ,HX,HY,pkf,T; ext = 20);

sn = 1
loc = [2000 2000 50]
ot = [0.05]
mt = [1 1 1 0 0 0]
wavelet = Ricker(pkf,model2.medium.dt)
waveform = wavelet
sou = initsource(sn, loc, ot, mt, waveform, model2)


for it in 1:model2.medium.nT
    addsou!(model2.wf,sou,it)
    extrap!(model2)
    if mod(it,300) == 0
      plt.show(plt.imshow(model2.wf.vx[:,:,23]))
  end
end
