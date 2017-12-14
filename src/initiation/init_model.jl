#=== elastic medium 2d ===#
function initmodel{T1,T2,T3,T4<:Real}(
    pvel::Union{String,<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float32,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},
    svel::Union{Void,String,<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},
    rho::T1,
    DZ::Union{<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},
    HX::T2,
    pkf::T3,
    T::T4;
    ext = 10,
    iflag = 2,
    dz = nothing,
    dx = nothing,
    dt= nothing,
    nT= nothing)

    medium = model(pvel, svel, rho, DZ, HX, pkf, T, ext, iflag, dx, dz,dt,nT)
    nwf = initnelwf(medium.nDZ, medium.nHX, ext, iflag)
    FDC = fdc(4)
# = unsplit PML boundary = #
    wf = initelwf(medium.nDZ, medium.nHX, ext, iflag)
    pml = BD(nwf, medium)
    calpara = wfpara(medium)
    return elmod2d(medium, wf, nwf, FDC, pml, calpara)
end


#=== acoustic medium 2d ===#
function initmodel{T1,T2,T3,T4<:Real}(
    pvel::Union{String,<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},
    rho::T1,
    DZ::Union{<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},
    HX::T2,
    pkf::T3,
    T::T4;
    ext = 10,
    iflag = 2,
    dz = nothing,
    dx = nothing,
    dt= nothing,
    nT= nothing)

    medium = model(pvel,rho, DZ, HX, pkf, T, ext, iflag, dx, dz,dt,nT)
    nwf = initnacwf(medium.nDZ, medium.nHX, ext, iflag)
    FDC = fdc(4)
# = unsplit PML boundary = #
    wf = initacwf(medium.nDZ, medium.nHX, ext, iflag)
    pml = BD(nwf, medium)
    calpara = wfpara(medium)
    return acmod2d(medium, wf, nwf, FDC, pml, calpara)
end






############ elastic medium 3d ############
function initmodel{T1,T2,T3,T4,T5<:Real}(
    pvel::Union{String,<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},
    svel::Union{Void,String,<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},
    rho::T1,
    DZ::Union{<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},
    HX::T2,
    HY::T3,
    pkf::T4,
    T::T5;
    ext = 10,
    iflag = 2,
    dz = nothing,
    dx = nothing,
    dy = nothing,
    dt= nothing,
    nT= nothing)

    medium = model(pvel, svel, rho, DZ, HX, HY, pkf, T, ext, iflag, dx, dy, dz, dt,nT)
    nwf = initnelwf(medium.nDZ, medium.nHX, medium.nHY, ext, iflag)
    FDC = fdc(4)
# = unsplit PML boundary = #
    wf = initelwf(medium.nDZ, medium.nHX, medium.nHY, ext, iflag)
    pml = BD(nwf, medium)
    calpara = wfpara(medium)
    return elmod3d(medium, wf, nwf, FDC, pml, calpara)
end







############ acoustic medium 3d ############
function initmodel{T1,T2,T3,T4,T5<:Real}(
    pvel::Union{String,<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},
    rho::T1,
    DZ::Union{<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},
    HX::T2,
    HY::T3,
    pkf::T4,
    T::T5;
    ext = 10,
    iflag = 2,
    dz = nothing,
    dx = nothing,
    dy = nothing,
    dt= nothing,
    nT= nothing)

    medium = model(pvel, rho, DZ, HX, HY, pkf, T, ext, iflag, dx, dy, dz, dt,nT)
    nwf = initnacwf(medium.nDZ, medium.nHX, medium.nHY, ext, iflag)
    FDC = fdc(4)
# = unsplit PML boundary = #
    wf = initacwf(medium.nDZ, medium.nHX, medium.nHY, ext, iflag)
    pml = BD(nwf, medium)
    calpara = wfpara(medium)
    return acmod3d(medium, wf, nwf, FDC, pml, calpara)
end
