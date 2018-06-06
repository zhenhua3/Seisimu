
function run!(model::elmod3d)

 path="/home/lzh/Dropbox/Zhenhua/Ongoing/Seisimu/deps/builds/el3d_openmp.so"
  ccall((:el3d_openmp,path),
  Void,
  (Ptr{Cdouble}, Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
   Ptr{Cdouble}, Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
   Ptr{Cdouble}, Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
   Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
   Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
   Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
   Cint, Cint, Cint,
   Ptr{Cdouble}, Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble},
   Ptr{Cdouble}, Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble},
   Ptr{Cdouble}, Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble},
   Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
   Cdouble, Cdouble, Cdouble, Cdouble, Cint,
   Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
   model.wf.vx, model.nwf.BDnvx[2], model.nwf.BDnvx[3], model.nwf.BDnvx[1],
   model.pml.PVxBTxx, model.pml.PVxBTxy, model.pml.PVxBTxz,
   model.wf.vy, model.nwf.BDnvy[2], model.nwf.BDnvy[3], model.nwf.BDnvy[1],
   model.pml.PVyBTxy, model.pml.PVyBTyy, model.pml.PVyBTyz,
   model.wf.vz, model.nwf.BDnvz[2], model.nwf.BDnvz[3], model.nwf.BDnvz[1],
   model.pml.PVzBTxz, model.pml.PVzBTyz, model.pml.PVzBTzz,
   model.wf.txx, model.pml.PTxxBVx, model.pml.PTxxBVy, model.pml.PTxxBVz,
   model.wf.tyy, model.pml.PTyyBVx, model.pml.PTyyBVy, model.pml.PTyyBVz,
   model.wf.tzz, model.pml.PTzzBVx, model.pml.PTzzBVy, model.pml.PTzzBVz,
   model.nwf.BDntpp[2], model.nwf.BDntpp[3], model.nwf.BDntpp[1],
   model.wf.txy, model.nwf.BDntxy[2], model.nwf.BDntxy[3], model.nwf.BDntxy[1],
   model.pml.PTxyBVx, model.pml.PTxyBVy,
   model.wf.tyz, model.nwf.BDntyz[2], model.nwf.BDntyz[3], model.nwf.BDntyz[1],
   model.pml.PTyzBVy, model.pml.PTyzBVz,
   model.wf.txz, model.nwf.BDntxz[2], model.nwf.BDntxz[3], model.nwf.BDntxz[1],
   model.pml.PTxzBVx, model.pml.PTxzBVz,
   model.medium.rho, model.medium.lambda, model.medium.mu, model.fdc,
   model.medium.dt, model.medium.dx, model.medium.dy, model.medium.dz, model.medium.ext,
   model.pml.bhalf, model.pml.ahalf, model.pml.bfull, model.pml.afull)
end

function run!(model::acmod3d)

 path="/home/lzh/Dropbox/Zhenhua/Ongoing/Seisimu/deps/builds/ac3d_openmp.so"

 ccall((:ac3d_openmp,path),
 Void,
 (Ptr{Cdouble}, Cint, Cint, Cint, Ptr{Cdouble},
  Ptr{Cdouble}, Cint, Cint, Cint, Ptr{Cdouble},
  Ptr{Cdouble}, Cint, Cint, Cint, Ptr{Cdouble},
  Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
  Cint, Cint, Cint,
  Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
  Cdouble, Cdouble, Cdouble, Cdouble, Cint,
  Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
  model.wf.vx, model.nwf.BDnvx[2], model.nwf.BDnvx[3], model.nwf.BDnvx[1],
  model.pml.PVxBTpp,
  model.wf.vy, model.nwf.BDnvy[2], model.nwf.BDnvy[3], model.nwf.BDnvy[1],
  model.pml.PVyBTpp,
  model.wf.vz, model.nwf.BDnvz[2], model.nwf.BDnvz[3], model.nwf.BDnvz[1],
  model.pml.PVzBTpp,
  model.wf.tpp, model.pml.PTppBVx, model.pml.PTppBVy, model.pml.PTppBVz,
  model.nwf.BDntpp[2], model.nwf.BDntpp[3], model.nwf.BDntpp[1],
  model.medium.rho, model.medium.lambda, model.fdc,
  model.medium.dt, model.medium.dx, model.medium.dy, model.medium.dz, model.medium.ext,
  model.pml.bhalf, model.pml.ahalf, model.pml.bfull, model.pml.afull)
end


function run!(model::acmod3d, snapshot_path::String; intvl=1, threadim=[16,16,4], AssignedStreamNum=6)

  acMemcpyGroup = acMemcpy(model,AssignedStreamNum)
  blockdim = [Int64(ceil(model.medium.BDnHX/threadim[1])),
        Int64(ceil(acMemcpyGroup.vx_PV_offset[1]/threadim[2])),
        Int64(ceil(model.medium.BDnDZ/threadim[3]))]
  RowNumber = TotalHToDDataSizeInBytes(model.nwf)
  ColNumber = Int64(ceil(model.medium.nT/intvl))
  dim = RowNumber*ColNumber



  path="/home/lzh/Dropbox/Zhenhua/Ongoing/Seisimu/deps/src/cuda/ac3d_cuda.so"

  ccall((:ac3d_cuda,path),
  Void,
  (Ptr{Cdouble}, Cint, Cint, Cint, Ptr{Cdouble},
   Ptr{Cdouble}, Cint, Cint, Cint, Ptr{Cdouble},
   Ptr{Cdouble}, Cint, Cint, Cint, Ptr{Cdouble},
   Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
   Cint, Cint, Cint,
   Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
   Cint, Cdouble, Cdouble, Cdouble, Cdouble, Cint,
   Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
   Ptr{UInt8}, Clonglong, Cint, Ptr{Clonglong}, Ptr{Clonglong},
   Cint, Cint, Cint,
   Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
   Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
   Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
   Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
   model.wf.vx, model.nwf.BDnvx[2], model.nwf.BDnvx[3], model.nwf.BDnvx[1],
   model.pml.PVxBTpp,
   model.wf.vy, model.nwf.BDnvy[2], model.nwf.BDnvy[3], model.nwf.BDnvy[1],
   model.pml.PVyBTpp,
   model.wf.vz, model.nwf.BDnvz[2], model.nwf.BDnvz[3], model.nwf.BDnvz[1],
   model.pml.PVzBTpp,
   model.wf.tpp, model.pml.PTppBVx, model.pml.PTppBVy, model.pml.PTppBVz,
   model.nwf.BDntpp[2], model.nwf.BDntpp[3], model.nwf.BDntpp[1],
   model.medium.rho, model.medium.lambda, model.fdc,
   model.medium.nT, model.medium.dt, model.medium.dx, model.medium.dy, model.medium.dz, model.medium.ext,
   model.pml.bhalf, model.pml.ahalf, model.pml.bfull, model.pml.afull,
   snapshot_path, dim, intvl, threadim, blockdim,
   acMemcpyGroup.AssignedStreamNum, acMemcpyGroup.TotalStreamNum, acMemcpyGroup.RegStreamDim,
   acMemcpyGroup.vx_PV_start,acMemcpyGroup.vy_PV_start,acMemcpyGroup.vz_PV_start,
   acMemcpyGroup.tpp_PV_start,acMemcpyGroup.rho_PV_start,
   acMemcpyGroup.vx_PV_offset,acMemcpyGroup.vy_PV_offset,acMemcpyGroup.vz_PV_offset,
   acMemcpyGroup.tpp_PV_offset,acMemcpyGroup.rho_PV_offset,
   acMemcpyGroup.vx_SS_start,acMemcpyGroup.vy_SS_start,acMemcpyGroup.vz_SS_start,
   acMemcpyGroup.tpp_SS_start,acMemcpyGroup.lambda_SS_start,
   acMemcpyGroup.vx_SS_offset,acMemcpyGroup.vy_SS_offset,acMemcpyGroup.vz_SS_offset,
   acMemcpyGroup.tpp_SS_offset,acMemcpyGroup.lambda_SS_offset)
end


function run!(model::elmod2d)

 path="/home/lzh/Dropbox/Zhenhua/Ongoing/Seisimu/deps/builds/el2d_openmp.so"

 ccall((:el2d_openmp,path),
 Void,
 (Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble},
  Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble},
  Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
  Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
  Cint, Cint,
  Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble},
  Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
  Cdouble, Cdouble, Cdouble, Cint,
  Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
  model.wf.vx, model.nwf.BDnvx[2], model.nwf.BDnvx[1],
  model.pml.PVxBTxx, model.pml.PVxBTxz,
  model.wf.vz, model.nwf.BDnvz[2], model.nwf.BDnvz[1],
  model.pml.PVzBTxz, model.pml.PVzBTzz,
  model.wf.txx, model.pml.PTxxBVx, model.pml.PTxxBVz,
  model.wf.tzz, model.pml.PTzzBVx, model.pml.PTzzBVz,
  model.nwf.BDntpp[2], model.nwf.BDntpp[1],
  model.wf.txz, model.nwf.BDntxz[2], model.nwf.BDntxz[1],
  model.pml.PTxzBVx, model.pml.PTxzBVz,
  model.medium.rho, model.medium.lambda, model.medium.mu, model.fdc,
  model.medium.dt, model.medium.dx, model.medium.dz, model.medium.ext,
  model.pml.bhalf, model.pml.ahalf, model.pml.bfull, model.pml.afull)
end

function run!(model::acmod2d)

 path="/home/lzh/Dropbox/Zhenhua/Ongoing/Seisimu/deps/builds/ac2d_openmp.so"

 ccall((:ac2d_openmp,path),
 Void,
 (Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble},
  Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble},
  Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
  Cint, Cint,
  Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
  Cdouble, Cdouble, Cdouble, Cint,
  Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
  model.wf.vx, model.nwf.BDnvx[2], model.nwf.BDnvx[1],
  model.pml.PVxBTpp,
  model.wf.vz, model.nwf.BDnvz[2], model.nwf.BDnvz[1],
  model.pml.PVzBTpp,
  model.wf.tpp, model.pml.PTppBVx, model.pml.PTppBVz,
  model.nwf.BDntpp[2], model.nwf.BDntpp[1],
  model.medium.rho, model.medium.lambda, model.fdc,
  model.medium.dt, model.medium.dx, model.medium.dz, model.medium.ext,
  model.pml.bhalf, model.pml.ahalf, model.pml.bfull, model.pml.afull)
end
