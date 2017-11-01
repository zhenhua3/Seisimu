type source{T1,T2<:Real}
    tp :: String  # source type
    loc :: Array{T1} # source's physical location in z and x axis
    nloc :: Array{Int64} # source's discrete location in z and x axis without PML boundary
    BDnloc :: Array{Int64} # source's discrete location in z and x axis with PML boundary. It is directly related to source type.
    ot :: T2 # source origin time
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

function SourType(tp::String, ext::Int64, snz::Int64, snx::Int64, nwf::nwf2d, iflag::Int64)
    # snx and snz : source discrete location in x and z direction

    BDnloc = Array{Int64,2}(4,5) #BDnloc[snz*snx,tp,Svx,Svz]

    if iflag == 1 #free surface

      if tp == "expl"
          BDnloc[:,1] = [(snz)+(snx+ext-1)*(nwf.BDntpp[1]);1;0;0] # Tpp
          BDnloc[:,2] = [(snz)+(snx+ext)*(nwf.BDnvx[1]);0;0;0] # Vx
          BDnloc[:,3] = [(snz+1)+(snx+ext)*(nwf.BDnvx[1]);0;0;0] # Vx
          BDnloc[:,4] = [(snz)+(snx+ext)*(nwf.BDnvz[1]);0;0;0] # Vz
          BDnloc[:,5] = [(snz)+(snx+1+ext)*(nwf.BDnvz[1]);0;0;0] # Vz
      elseif tp =="sfx"
          BDnloc[:,1] = [(snz)+(snx+ext-1)*(nwf.BDntpp[1]);0;0;0] # Tpp
          BDnloc[:,2] = [(snz)+(snx+ext)*(nwf.BDnvx[1]);0;1;0] # Vx
          BDnloc[:,3] = [(snz+1)+(snx+ext)*(nwf.BDnvx[1]);0;0;0] # Vx
          BDnloc[:,4] = [(snz)+(snx+ext)*(nwf.BDnvz[1]);0;0;0] # Vz
          BDnloc[:,5] = [(snz)+(snx+1+ext)*(nwf.BDnvz[1]);0;0;0] # Vz
      elseif tp =="sfz"
        BDnloc[:,1] = [(snz)+(snx+ext-1)*(nwf.BDntpp[1]);0;0;0] # Tpp
        BDnloc[:,2] = [(snz)+(snx+ext)*(nwf.BDnvx[1]);0;0;0] # Vx
        BDnloc[:,3] = [(snz+1)+(snx+ext)*(nwf.BDnvx[1]);0;0;0] # Vx
        BDnloc[:,4] = [(snz)+(snx+ext)*(nwf.BDnvz[1]);0;0;1] # Vz
        BDnloc[:,5] = [(snz)+(snx+1+ext)*(nwf.BDnvz[1]);0;0;0] # Vz
      elseif tp =="scx"
        BDnloc[:,1] = [(snz)+(snx+ext-1)*(nwf.BDntpp[1]);0;0;0] # Tpp
        BDnloc[:,2] = [(snz)+(snx+ext)*(nwf.BDnvx[1]);0;-1;0] # Vx
        BDnloc[:,3] = [(snz+1)+(snx+ext)*(nwf.BDnvx[1]);0;1;0] # Vx
        BDnloc[:,4] = [(snz)+(snx+ext)*(nwf.BDnvz[1]);0;0;0] # Vz
        BDnloc[:,5] = [(snz)+(snx+1+ext)*(nwf.BDnvz[1]);0;0;0] # Vz
      elseif tp =="scz"
        BDnloc[:,1] = [(snz)+(snx+ext-1)*(nwf.BDntpp[1]);0;0;0] # Tpp
        BDnloc[:,2] = [(snz)+(snx+ext)*(nwf.BDnvx[1]);0;0;0] # Vx
        BDnloc[:,3] = [(snz+1)+(snx+ext)*(nwf.BDnvx[1]);0;0;0] # Vx
        BDnloc[:,4] = [(snz)+(snx+ext)*(nwf.BDnvz[1]);0;0;1] # Vz
        BDnloc[:,5] = [(snz)+(snx+1+ext)*(nwf.BDnvz[1]);0;0;-1] # Vz
      elseif tp=="DC"
        BDnloc[:,1] = [(snz)+(snx+ext-1)*(nwf.BDntpp[1]);0;0;0] # Tpp
        BDnloc[:,2] = [(snz)+(snx+ext)*(nwf.BDnvx[1]);0;-1;0] # Vx
        BDnloc[:,3] = [(snz+1)+(snx+ext)*(nwf.BDnvx[1]);0;1;0] # Vx
        BDnloc[:,4] = [(snz)+(snx+ext)*(nwf.BDnvz[1]);0;0;1] # Vz
        BDnloc[:,5] = [(snz)+(snx+1+ext)*(nwf.BDnvz[1]);0;0;-1] # Vz
      end

    elseif iflag == 2 #unlimited medium

      if tp == "expl"
          BDnloc[:,1] = [(snz+ext)+(snx+ext-1)*(nwf.BDntpp[1]);1;0;0] # Tpp
          BDnloc[:,2] = [(snz+ext)+(snx+ext)*(nwf.BDnvx[1]);0;0;0] # Vx
          BDnloc[:,3] = [(snz+1+ext)+(snx+ext)*(nwf.BDnvx[1]);0;0;0] # Vx
          BDnloc[:,4] = [(snz+ext)+(snx+ext)*(nwf.BDnvz[1]);0;0;0] # Vz
          BDnloc[:,5] = [(snz+ext)+(snx+1+ext)*(nwf.BDnvz[1]);0;0;0] # Vz
      elseif tp =="sfx"
          BDnloc[:,1] = [(snz+ext)+(snx+ext-1)*(nwf.BDntpp[1]);0;0;0] # Tpp
          BDnloc[:,2] = [(snz+ext)+(snx+ext)*(nwf.BDnvx[1]);0;1;0] # Vx
          BDnloc[:,3] = [(snz+1+ext)+(snx+ext)*(nwf.BDnvx[1]);0;0;0] # Vx
          BDnloc[:,4] = [(snz+ext)+(snx+ext)*(nwf.BDnvz[1]);0;0;0] # Vz
          BDnloc[:,5] = [(snz+ext)+(snx+1+ext)*(nwf.BDnvz[1]);0;0;0] # Vz
      elseif tp =="sfz"
        BDnloc[:,1] = [(snz+ext)+(snx+ext-1)*(nwf.BDntpp[1]);0;0;0] # Tpp
        BDnloc[:,2] = [(snz+ext)+(snx+ext)*(nwf.BDnvx[1]);0;0;0] # Vx
        BDnloc[:,3] = [(snz+1+ext)+(snx+ext)*(nwf.BDnvx[1]);0;0;0] # Vx
        BDnloc[:,4] = [(snz+ext)+(snx+ext)*(nwf.BDnvz[1]);0;0;1] # Vz
        BDnloc[:,5] = [(snz+ext)+(snx+1+ext)*(nwf.BDnvz[1]);0;0;0] # Vz
      elseif tp =="scx"
        BDnloc[:,1] = [(snz+ext)+(snx+ext-1)*(nwf.BDntpp[1]);0;0;0] # Tpp
        BDnloc[:,2] = [(snz+ext)+(snx+ext)*(nwf.BDnvx[1]);0;-1;0] # Vx
        BDnloc[:,3] = [(snz+1+ext)+(snx+ext)*(nwf.BDnvx[1]);0;1;0] # Vx
        BDnloc[:,4] = [(snz+ext)+(snx+ext)*(nwf.BDnvz[1]);0;0;0] # Vz
        BDnloc[:,5] = [(snz+ext)+(snx+1+ext)*(nwf.BDnvz[1]);0;0;0] # Vz
      elseif tp =="scz"
        BDnloc[:,1] = [(snz+ext)+(snx+ext-1)*(nwf.BDntpp[1]);0;0;0] # Tpp
        BDnloc[:,2] = [(snz+ext)+(snx+ext)*(nwf.BDnvx[1]);0;0;0] # Vx
        BDnloc[:,3] = [(snz+1+ext)+(snx+ext)*(nwf.BDnvx[1]);0;0;0] # Vx
        BDnloc[:,4] = [(snz+ext)+(snx+ext)*(nwf.BDnvz[1]);0;0;-1] # Vz
        BDnloc[:,5] = [(snz+ext)+(snx+1+ext)*(nwf.BDnvz[1]);0;0;1] # Vz
      elseif tp=="DC"
        BDnloc[:,1] = [(snz+ext)+(snx+ext-1)*(nwf.BDntpp[1]);0;0;0] # Tpp
        BDnloc[:,2] = [(snz+ext)+(snx+ext)*(nwf.BDnvx[1]);0;-1;0] # Vx
        BDnloc[:,3] = [(snz+1+ext)+(snx+ext)*(nwf.BDnvx[1]);0;1;0] # Vx
        BDnloc[:,4] = [(snz+ext)+(snx+ext)*(nwf.BDnvz[1]);0;0;-1] # Vz
        BDnloc[:,5] = [(snz+ext)+(snx+1+ext)*(nwf.BDnvz[1]);0;0;1] # Vz
      end

    end
    return BDnloc
end

function CreateSource{T1,T2<:Real}(
    loc::Array{T1},
    ot::T2,
    tp::String,
    waveform::Array{Float64},
    medium::medium2d,
    nwf::nwf2d,
    iflag::Int64,
    ext::Int64)

    tpinput = tp
    nloc = [Int64(round(loc[1]/medium.dz)),  Int64(round(loc[2]/medium.dx))]
    BDnloc = SourType(tpinput, ext, nloc[1], nloc[2], nwf, iflag)
    otinput = ot
    not = Int64(round(otinput/medium.dt))
    SourceWaveform = CreateWaveform(1, not, medium.nT, waveform)

    Sou = source(tpinput, loc, nloc, BDnloc, otinput, not, SourceWaveform)
  end


function initsource{T1,T2<:Real}(
    sn::Int64,
    loc::Array{T1}, # location
    ot::Union{T2,Array{T2}}, # origin time
    tp::Union{String,Array{String,1}}, # type
    waveform::Array{Float64},
    model::Union{spmod2d,nspmod2d})

    sou = Array{source}(sn)
    for i in 1:sn
        sou[i] = CreateSource(
        loc[i,:],
        ot[i],
        tp[i],
        waveform[:,i],
        model.medium,
        model.nwf,
        model.medium.iflag,
        model.medium.ext)
    end
    return sou
end
