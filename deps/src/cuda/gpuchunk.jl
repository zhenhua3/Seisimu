# data_dim : Total length of the slowest dimension in memory, 2D in X direction, 3D in Y direction
# allowed_GPU_dim : Maximum dimension that can be stored in assigned GPU memory
# Num_stream : Stream number that initially defined according to GPU memory. Can be changed according to real data size.

# All streams compose a Group, whose size is the same as stream number
# precision : 8 for Float64; 4 for Float32.
# full : vx, vz, txx, tzz, txz
# half : vy, tyy, txy, tyz

function create_stream_chunk(data_dim::Int64, allowed_GPU_dim::Int64, Num_stream::Int64)
    Num_group = Int64(ceil(data_dim/allowed_GPU_dim))
    reggroup_dim = allowed_GPU_dim
    endgroup_dim = mod(data_dim,reggroup_dim)


    if data_dim <= allowed_GPU_dim
        Num_stream = Int64(ceil(data_dim/allowed_GPU_dim*Num_stream))
        reggroup_dim = data_dim
        endgroup_dim = 0

        chk_group_full = zeros(Int64,2);
        chk_group_half = zeros(Int64,2);

        chk_group_full[2] = reggroup_dim;
        chk_group_half[2] = reggroup_dim-1;

    elseif data_dim > allowed_GPU_dim
        if endgroup_dim == 1 && Num_group > 2;
            endgroup_dim = reggroup_dim +1;
            Num_group = Num_group -1;

            chk_group_full = zeros(Int64,Num_group+1);
            chk_group_half = zeros(Int64,Num_group+1);

            for i in 2:Num_group
                chk_group_full[i] = reggroup_dim;
                chk_group_half[i] = reggroup_dim;
            end
            chk_group_full[Num_group+1] = endgroup_dim;
            chk_group_half[Num_group+1] = endgroup_dim-1;

        elseif endgroup_dim == 1 && Num_group == 2;
            Num_group == 1;
            endgroup_dim = 0;
            reggroup_dim = reggroup_dim + 1;

            chk_group_full = zeros(Int64,2);
            chk_group_half = zeros(Int64,2);

            chk_group_full[2] = reggroup_dim;
            chk_group_half[2] = reggroup_dim-1;
        else
            chk_group_full = zeros(Int64,Num_group+1);
            chk_group_half = zeros(Int64,Num_group+1);

            for i in 2:Num_group
                chk_group_full[i] = reggroup_dim;
                chk_group_half[i] = reggroup_dim;
            end
            chk_group_full[Num_group+1] = endgroup_dim;
            chk_group_half[Num_group+1] = endgroup_dim-1;
        end
    end


    chk_stream_full = zeros(Int64,Num_stream+2,Num_group)
    chk_stream_half = zeros(Int64,Num_stream+2,Num_group)

    if endgroup_dim != 0
        Num_stream_endgroup = Int64(round(endgroup_dim/reggroup_dim*Num_stream))
        Dim_stream_regchunk_endgroup = Int64(ceil(endgroup_dim/Num_stream_endgroup))
        Dim_stream_endchunk_endgroup = mod(endgroup_dim,Dim_stream_regchunk_endgroup)
        Dim_stream_regchunk_reggroup = Int64(ceil(reggroup_dim/Num_stream))
        Dim_stream_endchunk_reggroup = mod(reggroup_dim,Dim_stream_regchunk_reggroup)

        if Dim_stream_endchunk_reggroup != 0 && Dim_stream_endchunk_endgroup != 0
            for i in 1:Num_group-1
                chk_stream_full[1,i] = Num_stream
                chk_stream_half[1,i] = Num_stream
                for j in 3:Num_stream+1
                    chk_stream_full[j,i] = Dim_stream_regchunk_reggroup;
                    chk_stream_half[j,i] = Dim_stream_regchunk_reggroup;
                end
                chk_stream_full[Num_stream+2,i] = Dim_stream_endchunk_reggroup;
                chk_stream_half[Num_stream+2,i] = Dim_stream_endchunk_reggroup;
            end
            i = Num_group
            if Dim_stream_endchunk_endgroup == 1
                Num_stream_endgroup = Num_stream_endgroup - 1 ;
                Dim_stream_endchunk_endgroup = Dim_stream_regchunk_endgroup + 1;
            end
            chk_stream_full[1,i] = Num_stream_endgroup
            chk_stream_half[1,i] = Num_stream_endgroup
            for j in 3:Num_stream_endgroup+1
                chk_stream_full[j,i] = Dim_stream_regchunk_endgroup;
                chk_stream_half[j,i] = Dim_stream_regchunk_endgroup;
            end
            chk_stream_full[Num_stream_endgroup+2,i] = Dim_stream_endchunk_endgroup;
            chk_stream_half[Num_stream_endgroup+2,i] = Dim_stream_endchunk_endgroup-1;

        elseif Dim_stream_endchunk_reggroup == 0 && Dim_stream_endchunk_endgroup != 0
            for i in 1:Num_group-1
                chk_stream_full[1,i] = Num_stream
                chk_stream_half[1,i] = Num_stream
                for j in 3:Num_stream+2
                    chk_stream_full[j,i] = Dim_stream_regchunk_reggroup;
                    chk_stream_half[j,i] = Dim_stream_regchunk_reggroup;
                end
            end
            i = Num_group
            if Dim_stream_endchunk_endgroup == 1
                Num_stream_endgroup = Num_stream_endgroup - 1 ;
                Dim_stream_endchunk_endgroup = Dim_stream_regchunk_endgroup + 1;
            end
            chk_stream_full[1,i] = Num_stream_endgroup
            chk_stream_half[1,i] = Num_stream_endgroup
            for j in 3:Num_stream_endgroup+1
                chk_stream_full[j,i] = Dim_stream_regchunk_endgroup;
                chk_stream_half[j,i] = Dim_stream_regchunk_endgroup;
            end
            chk_stream_full[Num_stream_endgroup+2,i] = Dim_stream_endchunk_endgroup;
            chk_stream_half[Num_stream_endgroup+2,i] = Dim_stream_endchunk_endgroup-1;

        elseif Dim_stream_endchunk_reggroup != 0 && Dim_stream_endchunk_endgroup == 0
            for i in 1:Num_group-1
                chk_stream_full[1,i] = Num_stream
                chk_stream_half[1,i] = Num_stream
                for j in 3:Num_stream+1
                    chk_stream_full[j,i] = Dim_stream_regchunk_reggroup;
                    chk_stream_half[j,i] = Dim_stream_regchunk_reggroup;
                end
                chk_stream_full[Num_stream+2,i] = Dim_stream_endchunk_reggroup;
                chk_stream_half[Num_stream+2,i] = Dim_stream_endchunk_reggroup;
            end
            i = Num_group
            chk_stream_full[1,i] = Num_stream_endgroup
            chk_stream_half[1,i] = Num_stream_endgroup
            for j in 3:Num_stream_endgroup+2
                chk_stream_full[j,i] = Dim_stream_regchunk_endgroup;
                chk_stream_half[j,i] = Dim_stream_regchunk_endgroup;
            end
            chk_stream_half[Num_stream_endgroup+2,Num_group] = Dim_stream_regchunk_endgroup-1;

        elseif Dim_stream_endchunk_reggroup == 0 && Dim_stream_endchunk_endgroup == 0

            for i in 1:Num_group-1
                chk_stream_full[1,i] = Num_stream
                chk_stream_half[1,i] = Num_stream
                for j in 3:Num_stream+2
                    chk_stream_full[j,i] = Dim_stream_regchunk_reggroup;
                    chk_stream_half[j,i] = Dim_stream_regchunk_reggroup;
                end
            end
            i = Num_group
            chk_stream_full[1,i] = Num_stream_endgroup
            chk_stream_half[1,i] = Num_stream_endgroup
            for j in 3:Num_stream_endgroup+2
                chk_stream_full[j,i] = Dim_stream_regchunk_endgroup;
                chk_stream_half[j,i] = Dim_stream_regchunk_endgroup;
            end
            chk_stream_half[Num_stream_endgroup+2,Num_group] = Dim_stream_regchunk_endgroup-1;
        end
    elseif endgroup_dim == 0
        Dim_stream_regchunk_reggroup = Int64(ceil(reggroup_dim/Num_stream))
        Dim_stream_endchunk_reggroup = mod(reggroup_dim,Dim_stream_regchunk_reggroup)

        if Dim_stream_endchunk_reggroup == 1 && Num_stream != 2
            Num_stream = Num_stream - 1 ;
            Dim_stream_endchunk_reggroup = Dim_stream_regchunk_reggroup + 1;
        elseif Dim_stream_endchunk_reggroup == 1 && Num_stream == 2
            Num_stream = 1
            Dim_stream_regchunk_reggroup = Dim_stream_regchunk_reggroup + 1;
            Dim_stream_endchunk_reggroup = 0;
        end

        if Dim_stream_endchunk_reggroup != 0
            for i in 1:Num_group
                chk_stream_full[1,i] = Num_stream
                chk_stream_half[1,i] = Num_stream
                for j in 3:Num_stream+1
                    chk_stream_full[j,i] = Dim_stream_regchunk_reggroup;
                    chk_stream_half[j,i] = Dim_stream_regchunk_reggroup;
                end
                chk_stream_full[Num_stream+2,i] = Dim_stream_endchunk_reggroup;
                chk_stream_half[Num_stream+2,i] = Dim_stream_endchunk_reggroup-1;
            end

        elseif Dim_stream_endchunk_reggroup == 0
            for i in 1:Num_group
                chk_stream_full[1,i] = Num_stream
                chk_stream_half[1,i] = Num_stream
                for j in 3:Num_stream+2
                    chk_stream_full[j,i] = Dim_stream_regchunk_reggroup;
                    chk_stream_half[j,i] = Dim_stream_regchunk_reggroup;
                end
            end
            chk_stream_half[Num_stream+2,Num_group] = chk_stream_half[Num_stream+2,Num_group] - 1;
        end
    end
    return chk_group_full, chk_group_half, chk_stream_full, chk_stream_half
end

#2D
function gpuchunk{T1,T2<:Real}(Medium_Type::String,BDnDZ::Int64,BDnHX::Int64,
    Data_Type::String, GPU_mem::T1, CPU_mem::T2, Num_stream::Int64)

    if Data_Type == "double"
        precision = 8
    elseif Data_Type == "single"
        precision = 4
    else error("Choose from double or single")
    end

    if Num_stream == 0 error("Stream number can not be 0.") end
    if Medium_Type == "elastic"
        data_size = 8*BDnHX*BDnDZ*precision/1024/1024/1024
        if data_size + GPU_mem > CPU_mem
            error("Model is over limit for the best performamce. Reduce the size of model.")
        end
        efct_GPU_mem = GPU_mem - 8*BDnDZ*precision/1024/1024/1024
        allowed_GPU_dim = Int64(floor(efct_GPU_mem * 1024 * 1024 * 1024/8/BDnDZ/precision))
        if allowed_GPU_dim < 3
            error("Model is over the limit of GPU memory. Reduce the size of model or increase GPU memory.")
        else
            if BDnHX <= Int64(ceil(0.3 * allowed_GPU_dim))
                error("Model is too small for GPU. Suggest using CPU instead.")
                chk_group_half, chk_group_full, chk_stream_half, chk_stream_full = create_stream_chunk(BDnHX-3,allowed_GPU_dim,Num_stream)
                return chk_group_full, chk_group_half, chk_stream_full, chk_stream_half
            else
                chk_group_half, chk_group_full, chk_stream_half, chk_stream_full = create_stream_chunk(BDnHX-3,allowed_GPU_dim,Num_stream)
            end
            return chk_group_full, chk_group_half, chk_stream_full, chk_stream_half
        end
    elseif Medium_Type == "acoustic"
        data_size = 6*BDnHX*BDnDZ*precision/1024/1024/1024
        if data_size + GPU_mem > CPU_mem
            error("Model is over limit for the best performamce. Reduce the size of model.")
        end
        efct_GPU_mem = GPU_mem - 4*BDnDZ*precision/1024/1024/1024
        allowed_GPU_dim = Int64(floor(efct_GPU_mem * 1024 * 1024 * 1024/6/BDnDZ/precision))
        if allowed_GPU_dim < 3
            error("Model is over the limit of GPU memory. Reduce the size of model or increase GPU memory.")
        else
            if BDnHX <= Int64(ceil(0.3 * allowed_GPU_dim))
                error("Model is too small for GPU. Suggest using CPU instead..")
                chk_group_half, chk_group_full, chk_stream_half, chk_stream_full = create_stream_chunk(BDnHX-3,allowed_GPU_dim,Num_stream)
                return chk_group_full, chk_group_half, chk_stream_full, chk_stream_half
            else
                chk_group_half, chk_group_full, chk_stream_half, chk_stream_full  = create_stream_chunk(BDnHX-3,allowed_GPU_dim,Num_stream)
                return chk_group_full, chk_group_half, chk_stream_full, chk_stream_half
            end
        end
    end
end

#3D
function gpuchunk{T1,T2<:Real}(Medium_Type::String,BDnDZ::Int64,BDnHX::Int64,BDnHY::Int64,
    Data_Type::String,GPU_mem::T1,CPU_mem::T2, Num_stream::Int64)

    if Data_Type == "double"
        precision = 8
    elseif Data_Type == "single"
        precision = 4
    else error("Choose from double or single")
    end

    if Num_stream == 0 error("Stream number can not be 0.") end
    if Medium_Type == "elastic"
        data_size = 12*BDnHY*BDnHX*BDnDZ*precision/1024/1024/1024
        if data_size + GPU_mem > CPU_mem
            error("Model is over limit for the best performamce. Reduce the size of model.")
        end
        efct_GPU_mem = GPU_mem - 20*BDnHX*BDnDZ*precision/1024/1024/1024
        allowed_GPU_dim = Int64(floor(efct_GPU_mem * 1024 * 1024 * 1024/12/BDnDZ/BDnHX/precision))
        if allowed_GPU_dim < 3
            error("Model is over the limit of GPU memory. Reduce the size of model or increase GPU memory.")
        else
            if BDnHY <= Int64(ceil(0.3 * allowed_GPU_dim))
                error("Model is too small for GPU. Suggest using CPU instead.")
                chk_group_half, chk_group_full, chk_stream_half, chk_stream_full = create_stream_chunk(BDnHY-3,allowed_GPU_dim,Num_stream)
                return chk_group_full, chk_group_half, chk_stream_full, chk_stream_half
            else
                chk_group_half, chk_group_full, chk_stream_half, chk_stream_full = create_stream_chunk(BDnHY-3,allowed_GPU_dim,Num_stream)
                return chk_group_full, chk_group_half, chk_stream_full, chk_stream_half
            end
        end
    elseif Medium_Type == "acoustic"
        data_size = 8*BDnHY*BDnHX*BDnDZ*precision/1024/1024/1024
        if data_size + GPU_mem > CPU_mem
            error("Model is over limit for the best performamce. Reduce the size of model.")
        end
        efct_GPU_mem = GPU_mem - 13*BDnHX*BDnDZ*precision/1024/1024/1024
        allowed_GPU_dim = Int64(floor(efct_GPU_mem * 1024 * 1024 * 1024/8/BDnHX/BDnDZ/precision))
        if allowed_GPU_dim < 3
            error("Model is over the limit of GPU memory. Reduce the size of model or increase GPU memory.")
        else
            if BDnHY <= Int64(ceil(0.3 * allowed_GPU_dim))
                error("Model is too small for GPU. Suggest using CPU instead.")
                chk_group_half, chk_group_full, chk_stream_half, chk_stream_full = create_stream_chunk(BDnHY-3,allowed_GPU_dim,Num_stream)
                return chk_group_full, chk_group_half, chk_stream_full, chk_stream_half
            else
                chk_group_half, chk_group_full, chk_stream_half, chk_stream_full = create_stream_chunk(BDnHY-3,allowed_GPU_dim,Num_stream)
                return chk_group_full, chk_group_half, chk_stream_full, chk_stream_half
            end
        end
    end
end
