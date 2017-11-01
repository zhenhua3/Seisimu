type spmod2d
    medium::medium2d
    wf::spwf2d
    nwf::nwf2d
    fd::spfdmtx2d
end
type nspmod2d
    medium::medium2d
    wf::nspwf2d
    nwf::nwf2d
    fd::nspfdmtx2d
    pml::BD2d
end

############ elastic medium ############
function initmodel{T3,T5,T6,T7<:Real}(
    pvel::Union{String,<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},
    svel::Union{Void,String,<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},
    rho::T3,
    DZ::Union{<:Real,Array{<:Real,1},StepRange{Int64,Int64},UnitRange{Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}},
    HX::T5,
    pkf::T6,
    T::T7;
    ext = 10,
    iflag = 2,
    dz = nothing,
    dx = nothing,
    dt= nothing,
    nT= nothing,
    mor = false)

  medium = model(pvel, svel, rho, DZ, HX, pkf, T, ext, iflag, dx, dz,dt,nT)
  nwf = initnwf(medium.nDZ, medium.nHX, ext, iflag)
  FDC = FDCoeff(4)
# =================== unsplit PML boundary ===================================#
  if mor == false
    wf = initwf(medium.nDZ, medium.nHX, ext, iflag, mor)
    fd = fdmtx(nwf, medium, FDC, ext)
    pml = BD(nwf, medium, ext, iflag)
    return nspmod2d(medium, wf, nwf, fd, pml)
# =================== split PML boundary ===================================#
  elseif mor == true
    wf = initwf(medium.NDep, medium.NHor, ext, iflag, mor)
    fd = fdmtx(nwf, medium, FDC, ext, iflag)
    return spmod2d(medium, wf, nwf, fd)
  end
end
