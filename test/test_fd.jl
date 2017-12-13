using seisimu

vp = [5000.0]
vs = [3500.0]
rho = 2.8
depth = 2000
Horizon = 2000
pkf = 30
T = 1
ext = 10
iflag = 2
dx = nothing
dz = nothing
dt = nothing
nT = nothing

model = initmodel(vp,vs,rho,depth,Horizon,pkf,T);

model.wf.txz = rand(Float64,model.nwf.BDntxz[1],model.nwf.BDntxz[2]);

model.wf.txx = rand(Float64,model.nwf.BDntpp[1],model.nwf.BDntpp[2]);

model.wf.tzz = rand(Float64,model.nwf.BDntpp[1],model.nwf.BDntpp[2]);

@time begin
for i in 1:200
extrap2d!(model);
end
end
