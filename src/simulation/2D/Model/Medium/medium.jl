type medium2d{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12<:Real}
    pvel::Array{T1,2}
    svel::Array{T2,2}
    rho::Array{T3,2}
    lambda::Array{T4,2}
    mu::Array{T5,2}
    dx::T6
    dz::T7
    dt::T8
    DZ::T9
    HX::T10
    nDZ::Int64
    nHX::Int64
    BDnDZ::Int64
    BDnHX::Int64
    T::T11 # Total Simulation Time
    nT::Int64 # Total discretized sampling time
    pkf::T12
    ext::Int64
    iflag::Int64
end

#=== add absorbing boundary ===#
function ModExpand(para::AbstractArray, ext::Int64, iflag::Int64)
            (nDZ,nHX) = size(para)
    if iflag == 1 #free surface
                temp_vertical2 = repmat(para[end:end,:],ext,1)
                para = vcat(para, temp_vertical2)
                temp_horizon1 = repmat(para[:,1],1,ext)
                temp_horizon2 = repmat(para[:,end],1,ext)
                para = hcat(temp_horizon2,para,temp_horizon1)
                BDnDZ = nDZ + ext
                BDnHX = nHX + 2*ext
        elseif iflag == 2  # unlimited medium
                temp_vertical1 = repmat(para[1:1,:],ext,1)
                temp_vertical2 = repmat(para[end:end,:],ext,1)
                para = vcat(temp_vertical1,para,temp_vertical2)
                temp_horizon1 = repmat(para[:,1],1,ext)
                temp_horizon2 = repmat(para[:,end],1,ext)
                para = hcat(temp_horizon1,para,temp_horizon2)
                BDnDZ = nDZ + 2*ext
                BDnHX = nHX + 2*ext
            end
        return para, BDnDZ, BDnHX
end







############ elastic layered medium ############
function model{T1,T2,T3,T4<:Real}(
    pvel::Union{String,<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},
    svel::Union{Void,String,<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},
    rho::T1,
    DZ::Union{<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},
    HX::T2,
    pkf::T3,
    T::T4,
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
        pvel = reshape(pvel, DZ,HX)
        svel = reshape(svel, DZ,HX)
        rho = zeros(DZ,HX) .+ rho
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
        Pvel = zeros(nDZ,nHX)
        Svel = zeros(nDZ,nHX)
        rho = zeros(nDZ,nHX) .+ rho
        # velocity interval and depth interval

        Pvel[1:nlayerdep[1],:] = pvel[1]
        Svel[1:nlayerdep[1],:] = svel[1]
        SumDep = 0
        for iL = 2 : length(pvel)
            SumDep = SumDep + nlayerdep[iL-1]
            Pvel[SumDep+1:SumDep+nlayerdep[iL],:] = pvel[iL]
            Svel[SumDep+1:SumDep+nlayerdep[iL],:] = svel[iL]
        end
        pvel = Pvel
        svel = Svel
    end

    pvel, BDnDZ, BDnHX = ModExpand(pvel, ext, iflag) # 1 : free surface; 2: unlimited surface
    svel, BDnDZ, BDnHX = ModExpand(svel, ext, iflag)
     rho, BDnDZ, BDnHX = ModExpand(rho, ext, iflag)
    lambda = pvel.^2.*rho - 2*svel.^2.*rho
    mu = svel.^2.*rho
    return medium2d(pvel,svel,rho,lambda,mu,dx,dz,dt,DZ,HX,nDZ,nHX,BDnDZ,BDnHX,T,nT,pkf,ext,iflag)
end
