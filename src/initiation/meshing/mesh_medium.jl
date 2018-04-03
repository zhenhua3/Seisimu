#=== Absorbing Boundary ===#
function ModExpand(para::AbstractArray, ext::Int64, iflag::Union{Void,Int64})
    if iflag == 1 #free surface
        temp_vertical2 = repmat(para[end:end,:],ext,1)
        para = vcat(para, temp_vertical2)
        temp_horizon1 = repmat(para[:,1],1,ext)
        temp_horizon2 = repmat(para[:,end],1,ext)
        para = hcat(temp_horizon2,para,temp_horizon1)
    elseif iflag == 2 # unlimited medium
        temp_vertical1 = repmat(para[1:1,:],ext,1)
        temp_vertical2 = repmat(para[end:end,:],ext,1)
        para = vcat(temp_vertical1,para,temp_vertical2)
        temp_horizon1 = repmat(para[:,1],1,ext)
        temp_horizon2 = repmat(para[:,end],1,ext)
        para = hcat(temp_horizon1,para,temp_horizon2)
    elseif iflag == nothing
        temp_horizon1 = repmat(para[:,1],1,ext)
        temp_horizon2 = repmat(para[:,end],1,ext)
        para = hcat(temp_horizon1,para,temp_horizon2)
    end
    return para
end






############ elastic medium 2d ############
function model{T1,T2,T3<:Real}(
    pvel::Union{String,<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64}}},
    svel::Union{Void,String,<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64}}},
    rho::Union{<:Real,String},
    DZ::Union{<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},
    HX::T1,
    pkf::T2,
    T::T3,
    ext::Int64,
    iflag::Int64,
    dx::Union{Void,<:Real},
    dz::Union{Void,<:Real},
    dt::Union{Void,<:Real},
    nT::Union{Void,Int64})


    if typeof(pvel) == String
        if !(typeof(DZ) == Int64 && typeof(HX) == Int64)
            error("Depth and Horzion have to be Integer grid number")
        elseif dx == nothing && dz == nothing
            error("Set up grid size dx and dz")
        elseif !(typeof(DZ) == Int64 && typeof(HX) == Int64)&&dx == nothing && dz == nothing
            error("1. Depth and Horzion have to be Integer grid number. 2.Set up grid size dx and dz.")
        end
        if pvel[end-2:end] == "bin"
            fid = open(pvel,"r")
            pvel = read(fid,Float64, DZ*HX,1)
            close(fid)
        elseif pvel[end-2:end] == "dat"
            fid = open(pvel,"r")
            pvel = readdlm(fid)
            close(fid)
        else error("P velocity file is not valid")
        end
        if svel == nothing
            svel = pvel./sqrt(3);
        elseif svel[end-2:end] == "bin"
            fid = open(svel,"r")
            svel = read(fid,Float64, DZ*HX,1)
            close(fid)
        elseif svel[end-2:end] == "dat"
            fid = open(svel,"r")
            svel = readdlm(fid)
            close(fid)
        else error("S velocity file is not valid")
        end
        max_pvel = maximum(pvel)
        min_pvel = minimum(pvel)
        max_svel = maximum(svel)
        min_svel = minimum(svel)

        if !(pkf < 2*min_svel/3/dx && pkf < 2*min_svel/3/dz)
            error("Peak frequency is too high for selected grid size")
        end
        if dt == nothing
            dt = 0.406 * dz / max_pvel
            nT = Int64(round(T/dt))
        elseif dt !=nothing
            if dt > 0.6*dz/max_pvel
                error("Simulation does not satisfy stablity criterion:
                maximum(pvel)*dt/dz < 0.6.")
            end
            if nT != nothing
                T= nT*dt
            elseif nT == nothing
                nT = Int64(round(T/dt))
            end
        end
        pvel = Float64.(reshape(pvel, DZ,HX))
        svel = Float64.(reshape(svel, DZ,HX))
        if typeof(rho) == String
            if rho[end-2:end] == "bin"
                fid = open(rho,"r")
                rho = read(fid,Float64, DZ*HX,1)
                close(fid)
            elseif rho[end-2:end] == "dat"
                fid = open(rho,"r")
                rho = readdlm(rho)
                close(fid)
            else error("Density file is not valid")
            end
            Rho = Float64.(reshape(rho, DZ,HX))
        elseif typeof(rho)<:Real
            Rho = Float64.(zeros(DZ,HX) .+ rho)
        else error("Give a density value or file path")
        end
        nDZ = DZ
        nHX = HX
        DZ = DZ*dz
        HX = HX*dx





    elseif typeof(pvel) != String
        if svel == nothing
            svel = pvel./sqrt(3);
        end
        if !(length(pvel)==length(svel)==length(DZ))
            error("Layer numbers of pvel, svel and depth should be the same. Please give a depth for each layer.")
        end
        max_pvel = maximum(pvel)
        min_pvel = minimum(pvel)
        max_svel = maximum(svel)
        min_svel = minimum(svel)
        # Spatial interval and time interval
        if !(max_svel > 0 && max_pvel > 0 && min_svel > 0 && max_pvel > 0)
            error("Medium velocities have to be larger than 0")
        else
            dz = min_svel / 20 / pkf
            dx = min_svel / 20 / pkf
        end
        if dt == nothing
            dt = 0.406 * dz / max_pvel
            nT = Int64(round(T/dt))
        elseif dt !=nothing
            if dt > 0.6*dz/max_pvel
                error("Simulation does not satisfy stablity criterion:
                maximum(P velocity)*dt/dz < 0.6.")
            end
            if nT != nothing
                T= nT*dt
            elseif nT == nothing
                nT = Int64(round(T/dt))
            end
        end
        nDZ = 0
        nlayerdep = zeros(Int64,length(DZ))
        if typeof(DZ) <: Array{<:Real,1} || typeof(DZ) <: Real
            for i = 1:length(DZ)
                nlayerdep[i] = Int64(round(DZ[i]/dz))
                nDZ = nDZ + nlayerdep[i]
            end
            DZ = sum(DZ)
        else
            nlayerdep[1] = Int64(round(DZ[1]/dz))
            nDZ = Int64(round(DZ[1]/dz))
            for i = 2:length(DZ)
                nlayerdep[i] = Int64(round((DZ[i]-DZ[i-1])/dz))
                nDZ = nDZ + nlayerdep[i]
            end
            DZ = DZ[end]
        end

        nHX = Int64(round(HX/dx))
        Pvel = zeros(Float64,nDZ,nHX)
        Svel = zeros(Float64,nDZ,nHX)
        if typeof(rho) == String
            if rho[end-2:end] == "bin"
                fid = open(rho,"r")
                rho = read(fid,Float64, DZ*HX,1)
                close(fid)
            elseif rho[end-2:end] == "dat"
                fid = open(rho,"r")
                rho = readdlm(rho)
                close(fid)
            else error("Density file is not valid")
            end
            Rho = Float64.(reshape(rho, DZ,HX))
        elseif typeof(rho)<:Real
            Rho = Float64.(zeros(DZ,HX) .+ rho)
        else error("Give a density value or file path")
        end

        # velocity interval and depth interval

        Pvel[1:nlayerdep[1],:] = Float64.(pvel[1])
        Svel[1:nlayerdep[1],:] = Float64.(svel[1])
        SumDep = 0
        for iL = 2 : length(pvel)
            SumDep = SumDep + nlayerdep[iL-1]
            Pvel[SumDep+1:SumDep+nlayerdep[iL],:] = Float64.(pvel[iL])
            Svel[SumDep+1:SumDep+nlayerdep[iL],:] = Float64.(svel[iL])
        end
        Pvel = Float64.(Pvel)
        Svel = Float64.(Svel)
    end

    pvel = ModExpand(Pvel, ext, iflag) # 1 : free surface; 2: unlimited surface
    svel = ModExpand(Svel, ext, iflag)
    rho = ModExpand(Rho, ext, iflag)
    BDnDZ,BDnHX = size(pvel)
    lambda = Array{Float64,2}(BDnDZ,BDnHX)
    mu = Array{Float64,2}(BDnDZ,BDnHX)
    Rho = Array{Float64,2}(BDnDZ,BDnHX)
    for i in 1 : BDnDZ
        for j in 1 : BDnHX
            lambda[i,j] = pvel[i,j]^2*rho[i,j]  - 2*svel[i,j]^2*rho[i,j]
            mu[i,j] = svel[i,j]^2*rho[i,j]
            Rho[i,j] = rho[i,j]
        end
    end

    return elastic2d(pvel,svel,Rho,lambda,mu,dx,dz,dt,DZ,HX,nDZ,nHX,BDnDZ,BDnHX,T,nT,pkf,ext,iflag)
end









############ acoustic medium 2d ############
  function model{T1,T2,T3<:Real}(
    pvel::Union{String,<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},
    rho::Union{<:Real,String},
    DZ::Union{<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},
    HX::T1,
    pkf::T2,
    T::T3,
    ext::Int64,
    iflag::Int64,
    dx::Union{Void,<:Real},
    dz::Union{Void,<:Real},
    dt::Union{Void,<:Real},
    nT::Union{Void,Int64})


    if typeof(pvel) == String
        if !(typeof(DZ) == Int64 && typeof(HX) == Int64)
            error("Depth and Horzion have to be Integer grid number")
        elseif dx == nothing && dz == nothing
            error("Set up grid size dx and dz")
        elseif !(typeof(DZ) == Int64 && typeof(HX) == Int64)&&dx == nothing && dz == nothing
            error("1. Depth and Horzion have to be Integer grid number. 2.Set up grid size dx and dz.")
        end
        if pvel[end-2:end] == "bin"
            fid = open(pvel,"r")
            pvel = read(fid,Float64, DZ*HX,1)
            close(fid)
        elseif pvel[end-2:end] == "dat"
            fid = open(pvel,"r")
            pvel = readdlm(fid)
            close(fid)
        else error("P velocity file is not valid")
        end

        max_pvel = maximum(pvel)
        min_pvel = minimum(pvel)

        if !(pkf < 2*min_pvel/3/dx && pkf < 2*min_pvel/3/dz)
            error("Peak frequency is too high for selected grid size")
        end
        if dt == nothing
            dt = 0.406 * dz / max_pvel
            nT = Int64(round(T/dt))
        elseif dt !=nothing
            if dt > 0.6*dz/max_pvel
                error("Simulation does not satisfy stablity criterion:
                maximum(pvel)*dt/dz < 0.6.")
            end
            if nT != nothing
                T= nT*dt
            elseif nT == nothing
                nT = Int64(round(T/dt))
            end
        end
        pvel = Float64.(reshape(pvel, DZ,HX))
        if typeof(rho) == String
            if rho[end-2:end] == "bin"
                fid = open(rho,"r")
                rho = read(fid,Float64, DZ*HX,1)
                close(fid)
            elseif rho[end-2:end] == "dat"
                fid = open(rho,"r")
                rho = readdlm(rho)
                close(fid)
            else error("Density file is not valid")
            end
            Rho = Float64.(reshape(rho, DZ,HX))
        elseif typeof(rho)<:Real
            Rho = Float64.(zeros(DZ,HX) .+ rho)
        else error("Give a density value or file path")
        end

        nDZ = DZ
        nHX = HX
        DZ = DZ*dz
        HX = HX*dx





    elseif typeof(pvel) != String

        if !(length(pvel)==length(DZ))
            error("Layer numbers of pvel, svel and depth should be the same. Please give a depth for each layer.")
        end
        max_pvel = maximum(pvel)
        min_pvel = minimum(pvel)
        # Spatial interval and time interval
        if !(max_pvel > 0 && max_pvel > 0)
            error("Medium velocities have to be larger than 0")
        else
            dz = min_pvel / 20 / pkf
            dx = min_pvel / 20 / pkf
        end
        if dt == nothing
            dt = 0.406 * dz / max_pvel
            nT = Int64(round(T/dt))
        elseif dt !=nothing
            if dt > 0.6*dz/max_pvel
                error("Simulation does not satisfy stablity criterion:
                maximum(P velocity)*dt/dz < 0.6.")
            end
            if nT != nothing
                T= nT*dt
            elseif nT == nothing
                nT = Int64(round(T/dt))
            end
        end
        nDZ = 0
        nlayerdep = zeros(Int64,length(DZ))
        if typeof(DZ) <: Array{<:Real,1} || typeof(DZ) <: Real
            for i = 1:length(DZ)
                nlayerdep[i] = Int64(round(DZ[i]/dz))
                nDZ = nDZ + nlayerdep[i]
            end
            DZ = sum(DZ)
        else
            nlayerdep[1] = Int64(round(DZ[1]/dz))
            nDZ = Int64(round(DZ[1]/dz))
            for i = 2:length(DZ)
                nlayerdep[i] = Int64(round((DZ[i]-DZ[i-1])/dz))
                nDZ = nDZ + nlayerdep[i]
            end
            DZ = DZ[end]
        end

        nHX = Int64(round(HX/dx))
        Pvel = zeros(Float64,nDZ,nHX)
        # velocity interval and depth interval

        Pvel[1:nlayerdep[1],:] = Float64.(pvel[1])
        SumDep = 0
        for iL = 2 : length(pvel)
            SumDep = SumDep + nlayerdep[iL-1]
            Pvel[SumDep+1:SumDep+nlayerdep[iL],:] = Float64.(pvel[iL])
        end
        Pvel = Float64.(Pvel)
    end

    if typeof(rho) == String
        if rho[end-2:end] == "bin"
            fid = open(rho,"r")
            rho = read(fid,Float64, DZ*HX,1)
            close(fid)
        elseif rho[end-2:end] == "dat"
            fid = open(rho,"r")
            rho = readdlm(rho)
            close(fid)
        else error("Density file is not valid")
        end
        Rho = Float64.(reshape(rho, DZ,HX))
    elseif typeof(rho)<:Real
        Rho = Float64.(zeros(DZ,HX) .+ rho)
    else error("Give a density value or file path")
    end

    pvel = ModExpand(Pvel, ext, iflag) # 1 : free surface; 2: unlimited surface
    rho = ModExpand(Rho, ext, iflag)

    BDnDZ, BDnHX = size(pvel)
    lambda = Array{Float64,2}(BDnDZ,BDnHX)
    Rho = Array{Float64,2}(BDnDZ,BDnHX)
    for i in 1 : BDnDZ
        for j in 1 : BDnHX
            lambda[i,j] = pvel[i,j]^2*rho[i,j]
            Rho[i,j] = rho[i,j]
        end
    end
    return acoustic2d(pvel,Rho,lambda,dx,dz,dt,DZ,HX,nDZ,nHX,BDnDZ,BDnHX,T,nT,pkf,ext,iflag)
end






############ elastic medium 3d ############
 function model{T1,T2,T3,T4<:Real}(
    pvel::Union{String,<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},
    svel::Union{Void,String,<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},
    rho::Union{<:Real,String},
    DZ::Union{<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},
    HX::T1,
    HY::T2,
    pkf::T3,
    T::T4,
    ext::Int64,
    iflag::Int64,
    dx::Union{Void,<:Real},
    dy::Union{Void,<:Real},
    dz::Union{Void,<:Real},
    dt::Union{Void,<:Real},
    nT::Union{Void,Int64})


    if typeof(pvel) == String
        if !(typeof(DZ) == Int64 && typeof(HX) == Int64 && typeof(HY) == Int64)
            error("Depth and Horzion have to be Integer grid number")
        elseif dx == nothing && dz == nothing && dy == nothing
            error("Set up grid size dx, dy and dz")
        elseif !(typeof(DZ) == Int64 && typeof(HX) == Int64 && typeof(HY) == Int64) && dx == nothing && dz == nothing && dy == nothing
            error("1. Depth and Horzion have to be Integer grid number. 2.Set up grid size dx, dy and dz.")
        end
        if pvel[end-2:end] == "bin"
            fid = open(pvel,"r")
            pvel = read(fid,Float64, DZ*HX*HY,1)
            close(fid)
        elseif pvel[end-2:end] == "dat"
            fid = open(pvel,"r")
            pvel = readdlm(fid)
            close(fid)
        else error("P velocity file is not valid")
        end
        if svel == nothing
            svel = pvel./sqrt(3);
        elseif svel[end-2:end] == "bin"
            fid = open(svel,"r")
            svel = read(fid,Float64, DZ*HX*HY,1)
            close(fid)
        elseif svel[end-2:end] == "dat"
            fid = open(svel,"r")
            svel = readdlm(fid)
            close(fid)
        else error("S velocity file is not valid")
        end
        max_pvel = maximum(pvel)
        min_pvel = minimum(pvel)
        max_svel = maximum(svel)
        min_svel = minimum(svel)

        if !(pkf < 2*min_svel/3/dx && pkf < 2*min_svel/3/dz && pkf < 2*min_svel/3/dy)
            error("Peak frequency is too high for selected grid size")
        end
        if dt == nothing
            dt = 0.406 * dz / max_pvel
            nT = Int64(round(T/dt))
        elseif dt != nothing
            if dt > 0.6 * dz / max_pvel
                error("Simulation does not satisfy stablity criterion:
                maximum(pvel)*dt/dz < 0.6.")
            end
            if nT != nothing
                T= nT*dt
            elseif nT == nothing
                nT = Int64(round(T/dt))
            end
        end
        Pvel = reshape(pvel, DZ,HX,HY)
        Svel = reshape(svel, DZ,HX,HY)
        if typeof(rho) == String
            if rho[end-2:end] == "bin"
                fid = open(rho,"r")
                rho = read(fid,Float64, DZ*HX*HY,1)
                close(fid)
            elseif rho[end-2:end] == "dat"
                fid = open(rho,"r")
                rho = readdlm(rho)
                close(fid)
            else error("Density file is not valid")
            end
            Rho = Float64.(reshape(rho, DZ,HX,HY))
        elseif typeof(rho)<:Real
            Rho = Float64.(zeros(DZ,HX,HY) .+ rho)
        else error("Give a density value or file path")
        end
        nDZ = DZ
        nHX = HX
        nHY = HY
        DZ = DZ*dz
        HX = HX*dx
        HY = HY*dy





    elseif typeof(pvel) != String
        if svel == nothing
            svel = pvel/sqrt(3);
        end
        if !(length(pvel)==length(svel)==length(DZ))
            error("Layer numbers of pvel, svel and depth should be the same. Please give a depth for each layer.")
        end
        max_pvel = maximum(pvel)
        min_pvel = minimum(pvel)
        max_svel = maximum(svel)
        min_svel = minimum(svel)
        # Spatial interval and time interval
        if !(max_svel > 0 && max_pvel > 0 && min_svel > 0 && max_pvel > 0)
            error("Medium velocities have to be larger than 0")
        else
            dz = min_svel / 20 / pkf
            dx = min_svel / 20 / pkf
            dy = min_svel / 20 / pkf
        end
        if dt == nothing
            dt = 0.406 * dz / max_pvel
            nT = Int64(round(T/dt))
        elseif dt != nothing
            if dt > 0.6*dz/max_pvel
                error("Simulation does not satisfy stablity criterion:
                maximum(P velocity)*dt/dz < 0.6.")
            end
            if nT != nothing
                T= nT*dt
            elseif nT == nothing
                nT = Int64(round(T/dt))
            end
        end
        nDZ = 0
        nlayerdep = zeros(Int64,length(DZ))
        if typeof(DZ) <: Array{<:Real,1} || typeof(DZ) <: Real
            for i = 1:length(DZ)
                nlayerdep[i] = Int64(round(DZ[i]/dz))
                nDZ = nDZ + nlayerdep[i]
            end
            DZ = sum(DZ)
        else
            nlayerdep[1] = Int64(round(DZ[1]/dz))
            nDZ = Int64(round(DZ[1]/dz))
            for i = 2:length(DZ)
                nlayerdep[i] = Int64(round((DZ[i]-DZ[i-1])/dz))
                nDZ = nDZ + nlayerdep[i]
            end
            DZ = DZ[end]
        end

        nHX = Int64(round(HX/dx))
        nHY = Int64(round(HY/dy))
        Pvel = zeros(nDZ,nHX,nHY)
        Svel = zeros(nDZ,nHX,nHY)
        if typeof(rho) == String
            if rho[end-2:end] == "bin"
                fid = open(rho,"r")
                rho = read(fid,Float64, DZ*HX*HY,1)
                close(fid)
            elseif rho[end-2:end] == "dat"
                fid = open(rho,"r")
                rho = readdlm(rho)
                close(fid)
            else error("Density file is not valid")
            end
            Rho = Float64.(reshape(rho, DZ,HX,HY))
        elseif typeof(rho)<:Real
            Rho = Float64.(zeros(DZ,HX,HY) .+ rho)
        else error("Give a density value or file path")
        end
        # velocity interval and depth interval

        Pvel[1:nlayerdep[1],:,:] = pvel[1]
        Svel[1:nlayerdep[1],:,:] = svel[1]
        SumDep = 0
        for iL = 2 : length(pvel)
            SumDep = SumDep + nlayerdep[iL-1]
            Pvel[SumDep+1:SumDep+nlayerdep[iL],:,:] = pvel[iL]
            Svel[SumDep+1:SumDep+nlayerdep[iL],:,:] = svel[iL]
        end
    end
    if iflag == 1 # free surface
        BDnDZ = nDZ + ext
    elseif iflag == 2
        BDnDZ = nDZ + 2*ext
    end
    BDnHX = nHX + 2*ext
    BDnHY = nHY + 2*ext
    pvel = zeros(BDnDZ,BDnHX,BDnHY)
    svel = zeros(BDnDZ,BDnHX,BDnHY)
     rho = zeros(BDnDZ,BDnHX,BDnHY)
    # 1 : free surface; 2: unlimited surface

    for i in 1:nHY
        pvel[:,:,ext+i] = ModExpand(Pvel[:,:,i], ext, iflag)
        svel[:,:,ext+i] = ModExpand(Svel[:,:,i], ext, iflag)
         rho[:,:,ext+i] = ModExpand( Rho[:,:,i], ext, iflag)
    end
    for i in 1:BDnDZ
        pvel[i,:,:] = ModExpand(pvel[i,:,ext+1:ext+nHY], ext, nothing)
        svel[i,:,:] = ModExpand(svel[i,:,ext+1:ext+nHY], ext, nothing)
         rho[i,:,:] = ModExpand( rho[i,:,ext+1:ext+nHY], ext, nothing)
    end
    # dz= Float64(dz)
    # dx= Float64(dx)
    # dy= Float64(dy)
    # dt= Float64(dt)
    lambda = Array{Float64,3}(BDnDZ,BDnHX,BDnHY)
    mu = Array{Float64,3}(BDnDZ,BDnHX,BDnHY)
    Rho = Array{Float64,3}(BDnDZ,BDnHX,BDnHY)
    for i in 1 : BDnDZ
        for j in 1 : BDnHX
            for k in 1 : BDnHY
                lambda[i,j,k] = pvel[i,j,k]^2*rho[i,j,k] - 2*svel[i,j,k]^2*rho[i,j,k]
                mu[i,j,k] = svel[i,j,k]^2*rho[i,j,k]
                Rho[i,j,k] = rho[i,j,k]
            end
        end
    end
    return elastic3d(pvel,svel,Rho,lambda,mu,dx,dy,dz,dt,
    DZ,HX,HY,nDZ,nHX,nHY,BDnDZ,BDnHX,BDnHY,T,nT,pkf,ext,iflag)
end









############ acoustic medium 3d ############
 function model{T1,T2,T3,T4<:Real}(
    pvel::Union{String,<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},
    rho::Union{<:Real,String},
    DZ::Union{<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},
    HX::T1,
    HY::T2,
    pkf::T3,
    T::T4,
    ext::Int64,
    iflag::Int64,
    dx::Union{Void,<:Real},
    dy::Union{Void,<:Real},
    dz::Union{Void,<:Real},
    dt::Union{Void,<:Real},
    nT::Union{Void,Int64})


    if typeof(pvel) == String
        if !(typeof(DZ) == Int64 && typeof(HX) == Int64 && typeof(HY) == Int64)
            error("Depth and Horzion have to be Integer grid number")
        elseif dx == nothing && dz == nothing && dy == nothing
            error("Set up grid size dx and dz")
        elseif !(typeof(DZ) == Int64 && typeof(HX) == Int64 && typeof(HY) == Int64) && dx == nothing && dz == nothing && dx == nothing
            error("1. Depth and Horzion have to be Integer grid number. 2.Set up grid size dx and dz.")
        end
        if pvel[end-2:end] == "bin"
            fid = open(pvel,"r")
            pvel = read(fid,Float64, DZ*HX*HY,1)
            close(fid)
        elseif pvel[end-2:end] == "dat"
            fid = open(pvel,"r")
            pvel = readdlm(fid)
            close(fid)
        else error("P velocity file is not valid")
        end

        max_pvel = maximum(pvel)
        min_pvel = minimum(pvel)

        if !(pkf < 2*min_pvel/3/dx && pkf < 2*min_pvel/3/dz && pkf < 2*min_pvel/3/dy)
            error("Peak frequency is too high for selected grid size")
        end
        if dt == nothing
            dt = 0.406 * dz / max_pvel
            nT = Int64(round(T/dt))
        elseif dt !=nothing
            if dt > 0.6*dz/max_pvel
                error("Simulation does not satisfy stablity criterion:
                maximum(pvel)*dt/dz < 0.6.")
            end
            if nT != nothing
                T= nT*dt
            elseif nT == nothing
                nT = Int64(round(T/dt))
            end
        end
        Pvel = reshape(pvel, DZ,HX,HY)
        if typeof(rho) == String
            if rho[end-2:end] == "bin"
                fid = open(rho,"r")
                rho = read(fid,Float64, DZ*HX*HY,1)
                close(fid)
            elseif rho[end-2:end] == "dat"
                fid = open(rho,"r")
                rho = readdlm(rho)
                close(fid)
            else error("Density file is not valid")
            end
            Rho = Float64.(reshape(rho, DZ,HX,HY))
        elseif typeof(rho)<:Real
            Rho = Float64.(zeros(DZ,HX,HY) .+ rho)
        else error("Give a density value or file path")
        end
        nDZ = DZ
        nHX = HX
        nHY = HY
        DZ = DZ*dz
        HX = HX*dx
        HY = HY*dy





    elseif typeof(pvel) != String

        if !(length(pvel)==length(DZ))
            error("Layer numbers of pvel, svel and depth should be the same. Please give a depth for each layer.")
        end
        max_pvel = maximum(pvel)
        min_pvel = minimum(pvel)
        # Spatial interval and time interval
        if !(max_pvel > 0 && max_pvel > 0)
            error("Medium velocities have to be larger than 0")
        else
            dz = min_pvel / 20 / pkf
            dx = min_pvel / 20 / pkf
            dy = min_pvel / 20 / pkf
        end
        if dt == nothing
            dt = 0.406 * dz / max_pvel
            nT = Int64(round(T/dt))
        elseif dt !=nothing
            if dt > 0.6*dz/max_pvel
                error("Simulation does not satisfy stablity criterion:
                maximum(P velocity)*dt/dz < 0.6.")
            end
            if nT != nothing
                T= nT*dt
            elseif nT == nothing
                nT = Int64(round(T/dt))
            end
        end
        nDZ = 0
        nlayerdep = zeros(Int64,length(DZ))
        if typeof(DZ) <: Array{<:Real,1} || typeof(DZ) <: Real
            for i = 1:length(DZ)
                nlayerdep[i] = Int64(round(DZ[i]/dz))
                nDZ = nDZ + nlayerdep[i]
            end
            DZ = sum(DZ)
        else
            nlayerdep[1] = Int64(round(DZ[1]/dz))
            nDZ = Int64(round(DZ[1]/dz))
            for i = 2:length(DZ)
                nlayerdep[i] = Int64(round((DZ[i]-DZ[i-1])/dz))
                nDZ = nDZ + nlayerdep[i]
            end
            DZ = DZ[end]
        end

        nHX = Int64(round(HX/dx))
        nHY = Int64(round(HY/dy))
        Pvel = zeros(nDZ,nHX,nHY)
        if typeof(rho) == String
            if rho[end-2:end] == "bin"
                fid = open(rho,"r")
                rho = read(fid,Float64, DZ*HX*HY,1)
                close(fid)
            elseif rho[end-2:end] == "dat"
                fid = open(rho,"r")
                rho = readdlm(rho)
                close(fid)
            else error("Density file is not valid")
            end
            Rho = Float64.(reshape(rho, DZ,HX,HY))
        elseif typeof(rho)<:Real
            Rho = Float64.(zeros(nDZ,nHX,nHY) .+ rho)
        else error("Give a density value or file path")
        end
        # velocity interval and depth interval
        for i in 1:nlayerdep[1]
            Pvel[i,:,:] = Pvel[i,:,:] .+ pvel[1]
        end
        SumDep = 0
        for iL = 2 : length(pvel)
            SumDep = SumDep + nlayerdep[iL-1]
            Pvel[SumDep+1:SumDep+nlayerdep[iL],:,:] = Pvel[SumDep+1:SumDep+nlayerdep[iL],:,:] .+ pvel[iL]
        end
    end

    # 1 : free surface; 2: unlimited surface
    if iflag == 1
        BDnDZ = nDZ + ext
    elseif iflag == 2
        BDnDZ = nDZ + 2*ext
    end
    BDnHX = nHX + 2*ext
    BDnHY = nHY + 2*ext
    pvel = zeros(BDnDZ,BDnHX,BDnHY)
     rho = zeros(BDnDZ,BDnHX,BDnHY)
    # 1 : free surface; 2: unlimited surface

    for i in 1:nHY
        pvel[:,:,ext+i] = ModExpand(Pvel[:,:,i], ext, iflag)
         rho[:,:,ext+i] = ModExpand( Rho[:,:,i], ext, iflag)
    end
    for i in 1:BDnDZ
        pvel[i,:,:] = ModExpand(pvel[i,:,ext+1:ext+nHY], ext, nothing)
         rho[i,:,:] = ModExpand( rho[i,:,ext+1:ext+nHY], ext, nothing)
    end

    lambda = Array{Float64,3}(BDnDZ,BDnHX,BDnHY)
    Rho = Array{Float64,3}(BDnDZ,BDnHX,BDnHY)
    for i in 1 : BDnDZ
        for j in 1 : BDnHX
            for k in 1 : BDnHY
                lambda[i,j,k] = pvel[i,j,k]^2*rho[i,j,k]
                Rho[i,j,k] = rho[i,j,k]
            end
        end
    end
    return acoustic3d(pvel,Rho,lambda,dx,dy,dz,dt,
    DZ,HX,HY,nDZ,nHX,nHY,BDnDZ,BDnHX,BDnHY,T,nT,pkf,ext,iflag)
end
