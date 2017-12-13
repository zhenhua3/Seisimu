type source{T1,T2, T3<:Real}
    mt :: Array{T1} # moment tensor
    loc :: Array{T2} # source's physical location in z and x axis
    nloc :: Array{Int64} # source's discrete location in z and x axis without PML boundary
    BDnloc :: Array{Int64} # source's discrete location in z and x axis with PML boundary. It is directly related to source type.
    ot :: T3 # source origin time
    not :: Int64 # source's discrete origin time
    waveform :: Array{Float64} # waveform for a single source
end

function CreateWaveform(sn::Int64, not::Int64, nT::Int64, waveform::Array{Float64})
  SourceWaveform = zeros(nT)
  twf = length(waveform)
  if not < nT && not > 0
    if not+twf-1 <= nT
      SourceWaveform[not:not+twf-1] = SourceWaveform[not:not+twf-1] + waveform
    elseif not+twf-1 > nT
      SourceWaveform[not:nT] = SourceWaveform[not:nT] + waveform[1:nT-not+1]
    end
  elseif not >= nT
    errormessage = string("check No.", string(sn), " source's starting time")
    error(errormessage)
  elseif not == 0
    if twf <= nT
      SourceWaveform[1:twf] = SourceWaveform[1:twf] + waveform
    elseif twf > nT
      SourceWaveform[1:nT] = SourceWaveform[1:nT] + waveform[1:nT]
    end
  end
  return SourceWaveform
end




function BDnlocation(ext::Int64, snz::Int64, snx::Int64, iflag::Int64)
    # snx and snz : source discrete location in x and z direction

    if iflag == 1 #free surface
        BDnloc = [snz, snx+ext]
    elseif iflag == 2 #unlimited medium
        BDnloc = [snz+ext, snx+ext]
    end
    return BDnloc
end


function BDnlocation(ext::Int64, snz::Int64, snx::Int64, sny::Int64, iflag::Int64)
    # snx and snz : source discrete location in x and z direction

    if iflag == 1 #free surface
        BDnloc = [snz, snx+ext, sny+ext]
    elseif iflag == 2 #unlimited medium
        BDnloc = [snz+ext, snx+ext, sny+ext]
    end
    return BDnloc
end



function CreateSource{T1,T2,T3<:Real}(
    loc :: Array{T1},
    ot :: T2,
    mt :: Array{T3,1}, # moment tensor (mxx,mzz,mxz)
    waveform :: Array{Float64},
    medium :: Union{elastic2d,acoustic2d},
    nwf :: Union{nelwf2d,nacwf2d})

    if length(mt) != 3
      error("Three elements in 2d moment tensor (mxx, mzz, mxz)")
    end

    nloc = [Int64(round(loc[1]/medium.dz)),  Int64(round(loc[2]/medium.dx))]
    BDnloc = BDnlocation(medium.ext, nloc[1], nloc[2], medium.iflag)
    otinput = ot
    not = Int64(round(otinput/medium.dt))
    SourceWaveform = CreateWaveform(1, not, medium.nT, waveform)

    Sou = source(mt, loc, nloc, BDnloc, otinput, not, SourceWaveform)
  end


  function CreateSource{T1,T2,T3<:Real}(
      loc::Array{T1},
      ot::T2,
      mt :: Array{T3,1}, # moment tensor (mxx,myy,mzz,mxy,mxz,myz)
      waveform::Array{Float64},
      medium::Union{elastic3d,acoustic3d},
      nwf::Union{nelwf3d,nacwf3d})

      if length(mt) != 6
        error("Six elements in 3d moment tensor (mxx, myy, mzz, mxy, mxz, myz)")
      end

      nloc = [Int64(round(loc[1]/medium.dz)),  Int64(round(loc[2]/medium.dx)),  Int64(round(loc[3]/medium.dy))]
      BDnloc = SourType(medium.ext, nloc[1], nloc[2], nloc[3], medium.iflag)
      otinput = ot
      not = Int64(round(otinput/medium.dt))
      SourceWaveform = CreateWaveform(1, not, medium.nT, waveform)

      Sou = source(mt, loc, nloc, BDnloc, otinput, not, SourceWaveform)
    end



#=== source ===#
function initsource{T1,T2, T3<:Real}(
    sn::Int64,
    loc::Array{T1}, # location
    ot::Union{T2,Array{T2}}, # origin time
    mt:: Array{T3}, # type
    waveform::Array{Float64},
    model::Union{elmod2d,acmod2d})

    sou = Array{source}(sn)
    for i in 1:sn
        sou[i] = CreateSource(
        loc[i,:],
        ot[i],
        mt[i,:],
        waveform[:,i],
        model.medium,
        model.nwf)
    end
    return sou
end


function initsource{T1,T2, T3<:Real}(
    sn::Int64,
    loc::Array{T1}, # location
    ot::Union{T2,Array{T2}}, # origin time
    mt:: Array{T3}, # type
    waveform::Array{Float64},
    model::Union{elmod3d,acmod3d})

    sou = Array{source}(sn)
    for i in 1:sn
        sou[i] = CreateSource(
        loc[i,:],
        ot[i],
        mt[i,:],
        waveform[:,i],
        model.medium,
        model.nwf)
    end
    return sou
end
