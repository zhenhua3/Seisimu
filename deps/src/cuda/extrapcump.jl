using seisimu

function extrap!{T1,T2<:Real}(model::elmod3d, threadim::Array{Int64}, CPU_mem::T1, GPU_mem::T2; num_stream=1)

    if length(threadim) != 3
        error("Thread dimension should be 3. Use [threadim1, threadim2, threadim3].")
    else
        chk_group_full, chk_group_half, chk_stream_full, chk_stream_half = gpuchunk("elastic",
        model.medium.BDnDZ, model.medium.BDnHX,
        model.medium.BDnHY,"double",GPU_mem, CPU_mem, num_stream)

        Num_group = length(chk_group_half[2:end])
        Max_group_dim = maximum(chk_group_half)
        Max_stream_dim = maximum(chk_stream_half[3:end,:])
        Max_num_stream = maximum(chk_stream_half[1,:])

        blockdim = [Int64(ceil(model.medium.BDnHX/threadim[1])),
        Int64(ceil(Max_stream_dim/threadim[2])),
        Int64(ceil(model.medium.BDnDZ/threadim[3]))]

        path="/home/lzh/Dropbox/Zhenhua/Ongoing/Seisimu/deps/src/cump/el3d_cump.so"
        ccall((:el3d_cump,path),
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
         Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
         Ptr{Cint}, Ptr{Cint}, Cint, Cint,
         Ptr{Cint}, Ptr{Cint}, Cint, Cint,
         Ptr{Cint}, Ptr{Cint}),
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
         model.pml.bhalf, model.pml.ahalf, model.pml.bfull, model.pml.afull,
         chk_group_full, chk_group_half, Num_group, Max_group_dim,
         chk_stream_full, chk_stream_half, Max_num_stream, Max_stream_dim,
         threadim, blockdim)
    end
end

function extrap!(model::acmod3d)

    path="/home/lzh/Dropbox/Zhenhua/Ongoing/Seisimu/deps/src/cuda/ac3d_cuda.so"


    ccall((:ac3d_cuda,path),
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
