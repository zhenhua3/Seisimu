function extrap2d!(model::elmod2d)

    path="/home/zhenhua3/Dropbox/Zhenhua/Ongoing/Seisimu/deps/builds/extrap2d_para"
    ccall((:el2dopmp,path),
    Void,
    (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
     Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
     Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
     Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
     Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
     Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Cint,
     Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
     model.wf.vx, model.wf.tmp_vx, model.nwf.nvx[1], model.nwf.nvx[2], model.calpara.rho_vx, model.pml.PVxBTxx, model.pml.PVxBTxz,
     model.wf.vz, model.wf.tmp_vz, model.nwf.nvz[1], model.nwf.nvz[2], model.calpara.rho_vz, model.pml.PVzBTzz, model.pml.PVzBTxz,
     model.wf.txx, model.wf.tzz, model.wf.tmp_tpp, model.nwf.ntpp[1], model.nwf.ntpp[2], model.calpara.lamu, model.medium.lambda, model.medium.mu,
     model.pml.PTxxBVx, model.pml.PTxxBVz, model.pml.PTzzBVx, model.pml.PTzzBVz,
     model.wf.txz, model.wf.tmp_txz, model.nwf.ntxz[1], model.nwf.ntxz[2], model.calpara.mu_txz, model.pml.PTxzBVx, model.pml.PTxzBVz,
     model.fdc, model.medium.dt, model.medium.dx, model.medium.dz, model.medium.ext,
     model.pml.bhalf, model.pml.ahalf, model.pml.bfull, model.pml.afull)

end

function extrap2d!(model::acmod2d)

    path="/home/zhenhua3/Dropbox/Zhenhua/Ongoing/Seisimu/deps/builds/extrap2d_para"
    ccall((:ac2dopmp,path),
    Void,
    (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble},
     Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble},
     Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble},
     Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
     Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Cint,
     Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
     model.wf.vx, model.wf.tmp_vx, model.nwf.nvx[1], model.nwf.nvx[2], model.calpara.rho_vx, model.pml.PVxBTxx,
     model.wf.vz, model.wf.tmp_vz, model.nwf.nvz[1], model.nwf.nvz[2], model.calpara.rho_vz, model.pml.PVzBTzz,
     model.wf.txx, model.wf.tzz, model.wf.tmp_tpp, model.nwf.ntpp[1], model.nwf.ntpp[2], model.medium.lambda,
     model.pml.PTxxBVx, model.pml.PTxxBVz, model.pml.PTzzBVx, model.pml.PTzzBVz,
     model.fdc, model.medium.dt, model.medium.dx, model.medium.dz, model.medium.ext,
     model.pml.bhalf, model.pml.ahalf, model.pml.bfull, model.pml.afull)

end
