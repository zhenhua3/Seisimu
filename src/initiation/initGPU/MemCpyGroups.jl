# data_dim : Total length of the slowest dimension in memory, 2D in X direction, 3D in Y direction
# allowed_GPU_dim : Maximum dimension that can be stored in assigned GPU memory
# Num_stream : Stream number that initially defined according to GPU memory. Can be changed according to real data size.

# All streams compose a Group, whose size is the same as stream number
# precision : 8 for Float64; 4 for Float32.
# full : vx, vz, txx, tzz, txz
# half : vy, tyy, txy, tyz

#3D
function MemCpyGroups(model::Union{elmod3d,acmod3d}, AssignedStreamNum::Int64)
    path="/home/lzh/Dropbox/Zhenhua/Ongoing/Seisimu/deps/src/cuda/ac3d_cuda.so"
    # available gpu memory in byte; 300 is extra space for data
    if typeof(model) == acmod3d
        AvailDeviceMem = ccall((:AvailDeviceMemInBytes,path),Cuint,())-6*3*AssignedStreamNum*model.medium.BDnDZ*model.medium.BDnHX*8
    elseif typeof(model) == elmod3d
        AvailDeviceMem = ccall((:AvailDeviceMemInBytes,path),Cuint,())-12*3*AssignedStreamNum*model.medium.BDnDZ*model.medium.BDnHX*8
    end

    TotalDataSize = TotalHToDDataSizeInBytes(model.nwf)
    if typeof(model) == elmod3d
        # 12 variables needed in device memory 8 bytes
        AllowedDeviceDim = Int64(floor(AvailDeviceMem/(12*model.medium.BDnDZ*model.medium.BDnHX*8)))
    elseif typeof(model) == acmod3d
        AllowedDeviceDim = Int64(floor(AvailDeviceMem/(6*model.medium.BDnDZ*model.medium.BDnHX*8)))
    end
    if AllowedDeviceDim < 5
        error("Model is over the limit of GPU memory. Reduce the size of model or increase GPU memory.")
    # elseif model.medium.BDnHY <= Int64(ceil(0.3 * AllowedDeviceDim))
    #     error("Model is too small for GPU. Suggest using CPU instead.")
    else
        if model.medium.BDnHY >= AllowedDeviceDim
            RegStreamDim = Int64(floor(AllowedDeviceDim/AssignedStreamNum))
        elseif model.medium.BDnHY < AllowedDeviceDim
            RegStreamDim = Int64(floor(model.medium.BDnHY/AssignedStreamNum))
        end
        if RegStreamDim == 0
            error("Assigned Stream Number is too large. Choose a smaller one.")
        else
            FinalStreamDim = mod(model.medium.BDnHY,RegStreamDim)
            if FinalStreamDim != 0
                TotalStreamNum = Int64((model.medium.BDnHY-FinalStreamDim)/RegStreamDim)+1
            elseif FinalStreamDim == 0
                TotalStreamNum = Int64(model.medium.BDnHY/RegStreamDim)
            end
            return TotalStreamNum, RegStreamDim, FinalStreamDim
        end
    end
end

function acMemcpy(model::acmod3d, AssignedStreamNum::Int64)
    TotalStreamNum, RegStreamDim, FinalStreamDim = MemCpyGroups(model, AssignedStreamNum)
    #SD Stream Dim
    # PV Particle Velocities
    # SS Stress
    vx_PV_offset = zeros(Int32,TotalStreamNum)
    vy_PV_offset = zeros(Int32,TotalStreamNum)
    vz_PV_offset = zeros(Int32,TotalStreamNum)
    tpp_PV_offset = zeros(Int32,TotalStreamNum)
    rho_PV_offset = zeros(Int32,TotalStreamNum)
    vx_PV_start = zeros(Int32,TotalStreamNum)
    vy_PV_start = zeros(Int32,TotalStreamNum)
    vz_PV_start = zeros(Int32,TotalStreamNum)
    tpp_PV_start = zeros(Int32,TotalStreamNum)
    rho_PV_start = zeros(Int32,TotalStreamNum)


    if TotalStreamNum == 1
        if RegStreamDim < 3 error("model is too small.") end
        vx_PV_offset[1] = vz_PV_offset[1] = tpp_PV_offset[1] = rho_PV_offset[1] = RegStreamDim
        vy_PV_offset[1] = RegStreamDim-2-1
        vx_PV_start[1] = vz_PV_start[1] = tpp_PV_start[1] = rho_PV_start[1] = 0
        vy_PV_start[1] = 1
    elseif TotalStreamNum == 2
        if FinalStreamDim == 0
            vx_PV_offset[1] = vz_PV_offset[1] = RegStreamDim
            vy_PV_offset[1] = RegStreamDim - 1
            tpp_PV_offset[1] = rho_PV_offset[1] = RegStreamDim + 2

            vx_PV_start[1] = vz_PV_start[1] = tpp_PV_start[1] = rho_PV_start[1] = 0
            vy_PV_start[1] = 1

            vx_PV_offset[2] = vz_PV_offset[2] = RegStreamDim
            vy_PV_offset[2] = RegStreamDim - 2
            tpp_PV_offset[2] = rho_PV_offset[2] = RegStreamDim + 1

            vx_PV_start[2] = vz_PV_start[2] = RegStreamDim
            vy_PV_start[2] = RegStreamDim
            tpp_PV_start[2] = rho_PV_start[2] = RegStreamDim - 1

        elseif FinalStreamDim <= 2 && FinalStreamDim > 0
            vx_PV_offset[1] = vz_PV_offset[1] = RegStreamDim
            vy_PV_offset[1] = RegStreamDim + FinalStreamDim - 3
            tpp_PV_offset[1] = rho_PV_offset[1] = RegStreamDim + FinalStreamDim

            vx_PV_start[1] = vz_PV_start[1] = tpp_PV_start[1] = rho_PV_start[1] = 0
            vy_PV_start[1] = 1

            vx_PV_offset[2] = vz_PV_offset[2] = FinalStreamDim
            vy_PV_offset[2] = 0
            tpp_PV_offset[2] = rho_PV_offset[2] = FinalStreamDim

            vx_PV_start[2] = vz_PV_start[2] = RegStreamDim
            vy_PV_start[2] = RegStreamDim
            tpp_PV_start[2] = rho_PV_start[2] = RegStreamDim
        elseif FinalStreamDim > 2
            vx_PV_offset[1] = vz_PV_offset[1] = RegStreamDim
            vy_PV_offset[1] = RegStreamDim - 1
            tpp_PV_offset[1] = rho_PV_offset[1] = RegStreamDim + 2

            vx_PV_start[1] = vz_PV_start[1] = tpp_PV_start[1] = rho_PV_start[1] = 0
            vy_PV_start[1] = 1

            vx_PV_offset[2] = vz_PV_offset[2] = FinalStreamDim
            vy_PV_offset[2] = FinalStreamDim - 1 - 1
            tpp_PV_offset[2] = rho_PV_offset[2] = FinalStreamDim + 1

            vx_PV_start[2] = vz_PV_start[2] = vy_PV_start[2] = RegStreamDim
            tpp_PV_start[2] = rho_PV_start[2] = RegStreamDim - 1
        end


    elseif TotalStreamNum >= 3
        vx_PV_offset[1] = vz_PV_offset[1] = RegStreamDim
        vy_PV_offset[1] = RegStreamDim - 1
        tpp_PV_offset[1] = rho_PV_offset[1] = RegStreamDim + 2

        vx_PV_start[1] = vz_PV_start[1] = tpp_PV_start[1] = rho_PV_start[1] = 0
        vy_PV_start[1] = 1

        for nstream in 2:TotalStreamNum-1
            vx_PV_offset[nstream] = vz_PV_offset[nstream] = vy_PV_offset[nstream] = RegStreamDim
            tpp_PV_offset[nstream] = rho_PV_offset[nstream] = RegStreamDim + 3

            vx_PV_start[nstream] = vz_PV_start[nstream] = vy_PV_start[nstream] = (nstream-1)*RegStreamDim
            tpp_PV_start[nstream] = rho_PV_start[nstream] = (nstream-1)*RegStreamDim - 1
        end

        if FinalStreamDim == 0

            vx_PV_offset[TotalStreamNum] = vz_PV_offset[TotalStreamNum] = RegStreamDim
            vy_PV_offset[TotalStreamNum] = RegStreamDim - 2
            tpp_PV_offset[TotalStreamNum] = rho_PV_offset[TotalStreamNum] = RegStreamDim + 1

            vx_PV_start[TotalStreamNum] = vz_PV_start[TotalStreamNum] = (TotalStreamNum-1)*RegStreamDim
            vy_PV_start[TotalStreamNum] = (TotalStreamNum-1)*RegStreamDim
            tpp_PV_start[TotalStreamNum] = rho_PV_start[TotalStreamNum] = (TotalStreamNum-1)*RegStreamDim - 1

        elseif FinalStreamDim <= 2 && FinalStreamDim > 0

            vx_PV_offset[TotalStreamNum-1] = vz_PV_offset[TotalStreamNum-1] = RegStreamDim
            vy_PV_offset[TotalStreamNum-1] = RegStreamDim + FinalStreamDim - 2
            tpp_PV_offset[TotalStreamNum-1] = rho_PV_offset[TotalStreamNum-1] = RegStreamDim + FinalStreamDim + 1

            vx_PV_start[TotalStreamNum-1] = vz_PV_start[TotalStreamNum-1] = (TotalStreamNum-2)*RegStreamDim
            vy_PV_start[TotalStreamNum-1] = (TotalStreamNum-2)*RegStreamDim
            tpp_PV_start[TotalStreamNum-1] = rho_PV_start[TotalStreamNum-1] = (TotalStreamNum-2)*RegStreamDim - 1

            vx_PV_offset[TotalStreamNum] = vz_PV_offset[TotalStreamNum] = FinalStreamDim
            vy_PV_offset[TotalStreamNum] = 0
            tpp_PV_offset[TotalStreamNum] = rho_PV_offset[TotalStreamNum] = FinalStreamDim

            vx_PV_start[TotalStreamNum] = vz_PV_start[TotalStreamNum] = (TotalStreamNum-1)*RegStreamDim
            vy_PV_start[TotalStreamNum] = (TotalStreamNum-1)*RegStreamDim
            tpp_PV_start[TotalStreamNum] = rho_PV_start[TotalStreamNum] = (TotalStreamNum-1)*RegStreamDim - 1

        elseif FinalStreamDim > 2

            vx_PV_offset[TotalStreamNum] = vz_PV_offset[TotalStreamNum] = FinalStreamDim
            vy_PV_offset[TotalStreamNum] = FinalStreamDim - 1 - 1
            tpp_PV_offset[TotalStreamNum] = rho_PV_offset[TotalStreamNum] = FinalStreamDim + 1

            vx_PV_start[TotalStreamNum] = vz_PV_start[TotalStreamNum] = (TotalStreamNum-1)*RegStreamDim
            vy_PV_start[TotalStreamNum] = (TotalStreamNum-1)*RegStreamDim
            tpp_PV_start[TotalStreamNum] = rho_PV_start[TotalStreamNum] = (TotalStreamNum-1)*RegStreamDim - 1

        end

    end

    vx_SS_offset = zeros(Int32,TotalStreamNum)
    vy_SS_offset = zeros(Int32,TotalStreamNum)
    vz_SS_offset = zeros(Int32,TotalStreamNum)
    tpp_SS_offset = zeros(Int32,TotalStreamNum)
    lambda_SS_offset = zeros(Int32,TotalStreamNum)
    vx_SS_start = zeros(Int32,TotalStreamNum)
    vy_SS_start = zeros(Int32,TotalStreamNum)
    vz_SS_start = zeros(Int32,TotalStreamNum)
    tpp_SS_start = zeros(Int32,TotalStreamNum)
    lambda_SS_start = zeros(Int32,TotalStreamNum)

    if TotalStreamNum == 1
        if RegStreamDim < 4 error("model is too small.") end
        tpp_SS_offset[1] = lambda_SS_offset[1] = RegStreamDim - 4
        vx_SS_offset[1] = vz_SS_offset[1] = RegStreamDim - 4
        vy_SS_offset[1] = RegStreamDim - 1

        vx_SS_start[1] = vz_SS_start[1] = 2
        tpp_SS_start[1] = lambda_SS_start[1] = 2
        vy_SS_start[1] = 0

    elseif TotalStreamNum == 2
        if FinalStreamDim == 0
            tpp_SS_offset[1] = lambda_SS_offset[1] = RegStreamDim - 2
            vx_SS_offset[1] = vz_SS_offset[1] = RegStreamDim - 2
            vy_SS_offset[1] = RegStreamDim + 1

            vx_SS_start[1] = vz_SS_start[1] = 2
            tpp_SS_start[1] = lambda_SS_start[1] = 2
            vy_SS_start[1] = 0

            tpp_SS_offset[2] = lambda_SS_offset[2] = RegStreamDim - 2
            vx_SS_offset[2] = vz_SS_offset[2] = RegStreamDim - 2
            vy_SS_offset[2] = RegStreamDim + 1

            vx_SS_start[2] = vz_SS_start[2] = RegStreamDim
            tpp_SS_start[2] = lambda_SS_start[2] = RegStreamDim
            vy_SS_start[2] = RegStreamDim - 2
        elseif FinalStreamDim <= 2
            tpp_SS_offset[1] = lambda_SS_offset[1] = RegStreamDim + FinalStreamDim - 4
            vx_SS_offset[1] = vz_SS_offset[1] = RegStreamDim + FinalStreamDim - 4
            vy_SS_offset[1] = RegStreamDim + FinalStreamDim - 1

            vx_SS_start[1] = vz_SS_start[1] = 2
            tpp_SS_start[1] = lambda_SS_start[1] = 2
            vy_SS_start[1] = 0

            tpp_SS_offset[2] = lambda_SS_offset[2] = 0
            vx_SS_offset[2] = vz_SS_offset[2] = 0
            vy_SS_offset[2] = 0

            vx_SS_start[2] = vz_SS_start[2] = RegStreamDim
            tpp_SS_start[2] = lambda_SS_start[2] = RegStreamDim
            vy_SS_start[2] = RegStreamDim - 2
        elseif FinalStreamDim > 2
            tpp_SS_offset[1] = lambda_SS_offset[1] = RegStreamDim - 2
            vx_SS_offset[1] = vz_SS_offset[1] = RegStreamDim - 2
            vy_SS_offset[1] = RegStreamDim + 1

            vx_SS_start[1] = vz_SS_start[1] = 2
            tpp_SS_start[1] = lambda_SS_start[1] = 2
            vy_SS_start[1] = 0

            tpp_SS_offset[2] = lambda_SS_offset[2] = FinalStreamDim - 2
            vx_SS_offset[2] = vz_SS_offset[2] = FinalStreamDim - 2
            vy_SS_offset[2] = FinalStreamDim + 1

            vx_SS_start[2] = vz_SS_start[2] = RegStreamDim
            tpp_SS_start[2] = lambda_SS_start[2] = RegStreamDim
            vy_SS_start[2] = RegStreamDim  - 2
        end

    elseif TotalStreamNum >= 3
        tpp_SS_offset[1] = lambda_SS_offset[1] = RegStreamDim - 2
        vx_SS_offset[1] = vz_SS_offset[1] = RegStreamDim - 2
        vy_SS_offset[1] = RegStreamDim + 1

        vx_SS_start[1] = vz_SS_start[1] = 2
        tpp_SS_start[1] = lambda_SS_start[1] = 2
        vy_SS_start[1] = 0

        for nstream in 2:TotalStreamNum-1
            tpp_SS_offset[nstream] = lambda_SS_offset[nstream] = RegStreamDim
            vx_SS_offset[nstream] = vz_SS_offset[nstream] = RegStreamDim
            vy_SS_offset[nstream] = RegStreamDim + 3

            vx_SS_start[nstream] = vz_SS_start[nstream] = (nstream-1)*RegStreamDim
            tpp_SS_start[nstream] = lambda_SS_start[nstream] = (nstream-1)*RegStreamDim
            vy_SS_start[nstream] = (nstream-1)*RegStreamDim - 2
        end

        if FinalStreamDim == 0
            tpp_SS_offset[TotalStreamNum] = lambda_SS_offset[TotalStreamNum] = RegStreamDim - 2
            vx_SS_offset[TotalStreamNum] = vz_SS_offset[TotalStreamNum] = RegStreamDim - 2
            vy_SS_offset[TotalStreamNum] = RegStreamDim + 1

            vx_SS_start[TotalStreamNum] = vz_SS_start[TotalStreamNum] = (TotalStreamNum-1)*RegStreamDim
            tpp_SS_start[TotalStreamNum] = lambda_SS_start[TotalStreamNum] = (TotalStreamNum-1)*RegStreamDim
            vy_SS_start[TotalStreamNum] = (TotalStreamNum-1)*RegStreamDim - 2

        elseif FinalStreamDim <= 2 && FinalStreamDim > 0
            tpp_SS_offset[TotalStreamNum-1] = lambda_SS_offset[TotalStreamNum-1] = RegStreamDim + FinalStreamDim - 2
            vx_SS_offset[TotalStreamNum-1] = vz_SS_offset[TotalStreamNum-1] = RegStreamDim + FinalStreamDim - 2
            vy_SS_offset[TotalStreamNum-1] = RegStreamDim + FinalStreamDim + 1

            vx_SS_start[TotalStreamNum-1] = vz_SS_start[TotalStreamNum-1] = (TotalStreamNum-2)*RegStreamDim
            tpp_SS_start[TotalStreamNum-1] = lambda_SS_start[TotalStreamNum-1] = (TotalStreamNum-2)*RegStreamDim
            vy_SS_start[TotalStreamNum-1] = (TotalStreamNum-2)*RegStreamDim - 2

            tpp_SS_offset[TotalStreamNum] = lambda_SS_offset[TotalStreamNum] = 0
            vx_SS_offset[TotalStreamNum] = vz_SS_offset[TotalStreamNum] = 0
            vy_SS_offset[TotalStreamNum] = 0

            vx_SS_start[TotalStreamNum] = vz_SS_start[TotalStreamNum] = (TotalStreamNum-1)*RegStreamDim
            tpp_SS_start[TotalStreamNum] = lambda_SS_start[TotalStreamNum] = (TotalStreamNum-1)*RegStreamDim
            vy_SS_start[TotalStreamNum] = (TotalStreamNum-1)*RegStreamDim - 2

        elseif FinalStreamDim > 2
            tpp_SS_offset[TotalStreamNum] = lambda_SS_offset[TotalStreamNum] = FinalStreamDim - 2
            vx_SS_offset[TotalStreamNum] = vz_SS_offset[TotalStreamNum] = FinalStreamDim - 2
            vy_SS_offset[TotalStreamNum] = FinalStreamDim + 1

            vx_SS_start[TotalStreamNum] = vz_SS_start[TotalStreamNum] = (TotalStreamNum-1)*RegStreamDim
            tpp_SS_start[TotalStreamNum] = lambda_SS_start[TotalStreamNum] = (TotalStreamNum-1)*RegStreamDim
            vy_SS_start[TotalStreamNum] = (TotalStreamNum-1)*RegStreamDim - 2
        end
    end

    return acstream(vx_PV_start,vy_PV_start,vz_PV_start,
                    tpp_PV_start,rho_PV_start,
                    vx_PV_offset,vy_PV_offset,vz_PV_offset,
                    tpp_PV_offset,rho_PV_offset,
                    vx_SS_start,vy_SS_start,vz_SS_start,
                    tpp_SS_start,lambda_SS_start,
                    vx_SS_offset,vy_SS_offset,vz_SS_offset,
                    tpp_SS_offset,lambda_SS_offset,
                    AssignedStreamNum, TotalStreamNum, RegStreamDim)
end
