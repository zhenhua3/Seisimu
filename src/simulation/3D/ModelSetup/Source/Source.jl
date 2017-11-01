type Source{T1,T2<:Real}
    stp :: String  # source type
    pos :: Array{T1} # source's physical location in z and x axis
    posn :: Array{Int64,1} # source's discrete location in z and x axis without PML boundary
    BDposn :: Array{Int64,2} # source's discrete location in z and x axis with PML boundary. It is directly related to source type.
    sot :: T2 # source origin time
    sotn :: Int64 # source's discrete origin time
    waveform :: Array{Float64} # waveform for a single source
end

type MultiSource{T1,T2<:Real}
  stp :: String # source type
  MultiSourPosition :: Array{T1,1} # physical position in m
  sot :: T2 # source origin time in s
  waveform :: Array{Float64,1}
end

function CreateWaveform(sn::Int64, sotn::Int64, Tn::Int64, waveform::Array{Float64})
  SourceWaveform = zeros(Tn)
  twf = length(waveform)
  if sotn < Tn && sotn > 0
    if sotn+twf-1 <= Tn
      SourceWaveform[sotn:sotn+twf-1] = SourceWaveform[sotn:sotn+twf-1] + waveform
    elseif sotn+twf-1 > Tn
      SourceWaveform[sotn:Tn] = SourceWaveform[sotn:Tn] + waveform[1:Tn-sotn+1]
    end
  elseif sotn >= Tn
    errormessage = string("check No.", string(sn), " source's starting time")
    error(errormessage)
  elseif sotn == 0
    if twf <= Tn
      SourceWaveform[1:twf] = SourceWaveform[1:twf] + waveform
    elseif twf > Tn
      SourceWaveform[1:Tn] = SourceWaveform[1:Tn] + waveform[1:Tn]
    end
  end
  return SourceWaveform
end

function MakeMultiSource{T1,T2<:Real}(sn::Int64, sourcetype::Array{String}, sourceposition::Array{T1,2}, sourceot::Array{T2,1}, sourwaveform::Array{Float64,2})
  sour = Array{MultiSource}(sn)
  for i = 1:sn
    sour[i] = MultiSource(sourcetype[i],sourceposition[:,i],sourceot[i],sourwaveform[:,i])
  end
  return sour
end

function SourType(stp::String, ext::Int64, snz::Int64, snx::Int64, WFSize::WF2DSize, iflag::Int64)
    # snx and snz : source discrete location in x and z direction

    BDposn = Array{Int64,2}(4,5) #BDposn[snz*snx,Stp,Svx,Svz]

    if iflag == 1 #free surface

      if stp == "expl"
          BDposn[:,1] = [(snz)+(snx+ext-1)*(WFSize.N_BDTpp[1]);1;0;0] # Tpp
          BDposn[:,2] = [(snz)+(snx+ext)*(WFSize.N_BDVx[1]);0;0;0] # Vx
          BDposn[:,3] = [(snz+1)+(snx+ext)*(WFSize.N_BDVx[1]);0;0;0] # Vx
          BDposn[:,4] = [(snz)+(snx+ext)*(WFSize.N_BDVz[1]);0;0;0] # Vz
          BDposn[:,5] = [(snz)+(snx+1+ext)*(WFSize.N_BDVz[1]);0;0;0] # Vz
      elseif stp =="sfx"
          BDposn[:,1] = [(snz)+(snx+ext-1)*(WFSize.N_BDTpp[1]);0;0;0] # Tpp
          BDposn[:,2] = [(snz)+(snx+ext)*(WFSize.N_BDVx[1]);0;1;0] # Vx
          BDposn[:,3] = [(snz+1)+(snx+ext)*(WFSize.N_BDVx[1]);0;0;0] # Vx
          BDposn[:,4] = [(snz)+(snx+ext)*(WFSize.N_BDVz[1]);0;0;0] # Vz
          BDposn[:,5] = [(snz)+(snx+1+ext)*(WFSize.N_BDVz[1]);0;0;0] # Vz
      elseif stp =="sfz"
        BDposn[:,1] = [(snz)+(snx+ext-1)*(WFSize.N_BDTpp[1]);0;0;0] # Tpp
        BDposn[:,2] = [(snz)+(snx+ext)*(WFSize.N_BDVx[1]);0;0;0] # Vx
        BDposn[:,3] = [(snz+1)+(snx+ext)*(WFSize.N_BDVx[1]);0;0;0] # Vx
        BDposn[:,4] = [(snz)+(snx+ext)*(WFSize.N_BDVz[1]);0;0;1] # Vz
        BDposn[:,5] = [(snz)+(snx+1+ext)*(WFSize.N_BDVz[1]);0;0;0] # Vz
      elseif stp =="scx"
        BDposn[:,1] = [(snz)+(snx+ext-1)*(WFSize.N_BDTpp[1]);0;0;0] # Tpp
        BDposn[:,2] = [(snz)+(snx+ext)*(WFSize.N_BDVx[1]);0;-1;0] # Vx
        BDposn[:,3] = [(snz+1)+(snx+ext)*(WFSize.N_BDVx[1]);0;1;0] # Vx
        BDposn[:,4] = [(snz)+(snx+ext)*(WFSize.N_BDVz[1]);0;0;0] # Vz
        BDposn[:,5] = [(snz)+(snx+1+ext)*(WFSize.N_BDVz[1]);0;0;0] # Vz
      elseif stp =="scz"
        BDposn[:,1] = [(snz)+(snx+ext-1)*(WFSize.N_BDTpp[1]);0;0;0] # Tpp
        BDposn[:,2] = [(snz)+(snx+ext)*(WFSize.N_BDVx[1]);0;0;0] # Vx
        BDposn[:,3] = [(snz+1)+(snx+ext)*(WFSize.N_BDVx[1]);0;0;0] # Vx
        BDposn[:,4] = [(snz)+(snx+ext)*(WFSize.N_BDVz[1]);0;0;1] # Vz
        BDposn[:,5] = [(snz)+(snx+1+ext)*(WFSize.N_BDVz[1]);0;0;-1] # Vz
      elseif stp=="DC"
        BDposn[:,1] = [(snz)+(snx+ext-1)*(WFSize.N_BDTpp[1]);0;0;0] # Tpp
        BDposn[:,2] = [(snz)+(snx+ext)*(WFSize.N_BDVx[1]);0;-1;0] # Vx
        BDposn[:,3] = [(snz+1)+(snx+ext)*(WFSize.N_BDVx[1]);0;1;0] # Vx
        BDposn[:,4] = [(snz)+(snx+ext)*(WFSize.N_BDVz[1]);0;0;1] # Vz
        BDposn[:,5] = [(snz)+(snx+1+ext)*(WFSize.N_BDVz[1]);0;0;-1] # Vz
      end

    elseif iflag == 2 #unlimited medium

      if stp == "expl"
          BDposn[:,1] = [(snz+ext)+(snx+ext-1)*(WFSize.N_BDTpp[1]);1;0;0] # Tpp
          BDposn[:,2] = [(snz+ext)+(snx+ext)*(WFSize.N_BDVx[1]);0;0;0] # Vx
          BDposn[:,3] = [(snz+1+ext)+(snx+ext)*(WFSize.N_BDVx[1]);0;0;0] # Vx
          BDposn[:,4] = [(snz+ext)+(snx+ext)*(WFSize.N_BDVz[1]);0;0;0] # Vz
          BDposn[:,5] = [(snz+ext)+(snx+1+ext)*(WFSize.N_BDVz[1]);0;0;0] # Vz
      elseif stp =="sfx"
          BDposn[:,1] = [(snz+ext)+(snx+ext-1)*(WFSize.N_BDTpp[1]);0;0;0] # Tpp
          BDposn[:,2] = [(snz+ext)+(snx+ext)*(WFSize.N_BDVx[1]);0;1;0] # Vx
          BDposn[:,3] = [(snz+1+ext)+(snx+ext)*(WFSize.N_BDVx[1]);0;0;0] # Vx
          BDposn[:,4] = [(snz+ext)+(snx+ext)*(WFSize.N_BDVz[1]);0;0;0] # Vz
          BDposn[:,5] = [(snz+ext)+(snx+1+ext)*(WFSize.N_BDVz[1]);0;0;0] # Vz
      elseif stp =="sfz"
        BDposn[:,1] = [(snz+ext)+(snx+ext-1)*(WFSize.N_BDTpp[1]);0;0;0] # Tpp
        BDposn[:,2] = [(snz+ext)+(snx+ext)*(WFSize.N_BDVx[1]);0;0;0] # Vx
        BDposn[:,3] = [(snz+1+ext)+(snx+ext)*(WFSize.N_BDVx[1]);0;0;0] # Vx
        BDposn[:,4] = [(snz+ext)+(snx+ext)*(WFSize.N_BDVz[1]);0;0;1] # Vz
        BDposn[:,5] = [(snz+ext)+(snx+1+ext)*(WFSize.N_BDVz[1]);0;0;0] # Vz
      elseif stp =="scx"
        BDposn[:,1] = [(snz+ext)+(snx+ext-1)*(WFSize.N_BDTpp[1]);0;0;0] # Tpp
        BDposn[:,2] = [(snz+ext)+(snx+ext)*(WFSize.N_BDVx[1]);0;-1;0] # Vx
        BDposn[:,3] = [(snz+1+ext)+(snx+ext)*(WFSize.N_BDVx[1]);0;1;0] # Vx
        BDposn[:,4] = [(snz+ext)+(snx+ext)*(WFSize.N_BDVz[1]);0;0;0] # Vz
        BDposn[:,5] = [(snz+ext)+(snx+1+ext)*(WFSize.N_BDVz[1]);0;0;0] # Vz
      elseif stp =="scz"
        BDposn[:,1] = [(snz+ext)+(snx+ext-1)*(WFSize.N_BDTpp[1]);0;0;0] # Tpp
        BDposn[:,2] = [(snz+ext)+(snx+ext)*(WFSize.N_BDVx[1]);0;0;0] # Vx
        BDposn[:,3] = [(snz+1+ext)+(snx+ext)*(WFSize.N_BDVx[1]);0;0;0] # Vx
        BDposn[:,4] = [(snz+ext)+(snx+ext)*(WFSize.N_BDVz[1]);0;0;-1] # Vz
        BDposn[:,5] = [(snz+ext)+(snx+1+ext)*(WFSize.N_BDVz[1]);0;0;1] # Vz
      elseif stp=="DC"
        BDposn[:,1] = [(snz+ext)+(snx+ext-1)*(WFSize.N_BDTpp[1]);0;0;0] # Tpp
        BDposn[:,2] = [(snz+ext)+(snx+ext)*(WFSize.N_BDVx[1]);0;-1;0] # Vx
        BDposn[:,3] = [(snz+1+ext)+(snx+ext)*(WFSize.N_BDVx[1]);0;1;0] # Vx
        BDposn[:,4] = [(snz+ext)+(snx+ext)*(WFSize.N_BDVz[1]);0;0;-1] # Vz
        BDposn[:,5] = [(snz+ext)+(snx+1+ext)*(WFSize.N_BDVz[1]);0;0;1] # Vz
      end

    end
    return BDposn
end

function source(sn::Int64, physour::Array{MultiSource}, ext::Int64,
    medium::Model, WFSize::WF2DSize, iflag::Int64)
# Multiple Source
# Input Parameters:
    # sn : number of sources
    # position : physical coordinate of sources, unit:m, is a matrix with size 2*sn, [z;x]
    # sot : source start time, unit: s
    # ext : extended PML boundary on one side
    # stp : source type, ASCIIString type Array,
    # medium :: Model type with all medium properties
    # waveform : Input wavelet for each source, size: nt*sn, nt: discrete wavelet sampling
    # T : total simulation time
    # iflag :: 1 for unlimited medium, 2 for free surface

    Sou = Array{Source}(sn)# compose source

    for i = 1:sn
        tpinput = physour[i].stp
        pos = physour[i].MultiSourPosition
        posn = [Int64(round(pos[1]/medium.dz)),  Int64(round(pos[2]/medium.dx))]
        BDposn = SourType(tpinput, ext, posn[1], posn[2], WFSize, iflag)
        otinput = physour[i].sot
        sotn = Int64(round(otinput/medium.dt))
        SourceWaveform = CreateWaveform(i, sotn, medium.Tn, physour[i].waveform)

        Sou[i] = Source(tpinput, pos, posn, BDposn, otinput, sotn, SourceWaveform)
    end

    return Sou
end

function source(sn::Int64, physour::Array{MultiSource}, ext::Int64, Tn::Int64,
    medium::Model, WFSize::WF2DSize, iflag::Int64)
# Multiple Source
# Input Parameters:
    # sn : number of sources
    # position : physical coordinate of sources, unit:m, is a matrix with size 2*sn, [z;x]
    # sot : source start time, unit: s
    # ext : extended PML boundary on one side
    # stp : source type, ASCIIString type Array,
    # medium :: Model type with all medium properties
    # waveform : Input wavelet for each source, size: nt*sn, nt: discrete wavelet sampling
    # T : total simulation time
    # iflag :: 1 for unlimited medium, 2 for free surface

    Sou = Array{Source}(sn)# compose source

    for i = 1:sn
        tpinput = physour[i].stp
        pos = physour[i].MultiSourPosition
        posn = [Int64(round(pos[1]/medium.dz)),  Int64(round(pos[2]/medium.dx))]
        BDposn = SourType(tpinput, ext, posn[1], posn[2], WFSize, iflag)
        otinput = physour[i].sot
        sotn = Int64(round(otinput/medium.dt))
        SourceWaveform = CreateWaveform(i, sotn, Tn, physour[i].waveform)

        Sou[i] = Source(tpinput, pos, posn, BDposn, otinput, sotn, SourceWaveform)
    end

    return Sou
end

function source{T1,T2<:Real}(position::Array{T1}, sot::T2,
    ext::Int64, stp::String,
    medium::Model, waveform::Array{Float64}, WFSize::WF2DSize, iflag::Int64)
    # Single source
    # Input Parameters:
        # position : physical coordinate of sources, unit:m, is a matrix with size 2*sn, [z;x]
        # sot : source start time, unit: s
        # ext : extended PML boundary on one side
        # stp : source type, ASCIIString type Array,
        # medium :: Model type with all medium properties
        # waveform : Input wavelet for each source, size: nt*sn, nt: discrete wavelet sampling
        # T : total simulation time
        # iflag :: 1 for unlimited medium, 2 for free surface

    tpinput = stp
    posn = [Int64(round(position[1]/medium.dz)),  Int64(round(position[2]/medium.dx))]
    BDposn = SourType(tpinput, ext, posn[1], posn[2], WFSize, iflag)
    otinput = sot
    sotn = Int64(round(otinput/medium.dt))
    SourceWaveform = CreateWaveform(1, sotn, medium.Tn, waveform)

    Sou = Source(tpinput, position, posn, BDposn, otinput, sotn, SourceWaveform)
  end

  function source{T1,T2<:Real}(position::Array{T1}, sot::T2,
      ext::Int64, Tn::Int64, stp::String,
      medium::Model, waveform::Array{Float64}, WFSize::WF2DSize, iflag::Int64)
      # Single source
      # Input Parameters:
          # position : physical coordinate of sources, unit:m, is a matrix with size 2*sn, [z;x]
          # sot : source start time, unit: s
          # ext : extended PML boundary on one side
          # stp : source type, ASCIIString type Array,
          # medium :: Model type with all medium properties
          # waveform : Input wavelet for each source, size: nt*sn, nt: discrete wavelet sampling
          # Tn : total simulation time
          # iflag :: 1 for unlimited medium, 2 for free surface

      tpinput = stp
      posn = [Int64(round(position[1]/medium.dz)),  Int64(round(position[2]/medium.dx))]
      BDposn = SourType(tpinput, ext, posn[1], posn[2], WFSize, iflag)
      otinput = sot
      sotn = Int64(round(otinput/medium.dt))
      SourceWaveform = CreateWaveform(1, sotn, Tn, waveform)

      Sou = Source(tpinput, position, posn, BDposn, otinput, sotn, SourceWaveform)
    end

function InitSou{T1,T2<:Real}(sn::Int64, position::Array{T1}, origintime::Union{T2,Array{T2,1}},
  sourcetype::Union{String,Array{String,1}}, waveform::Array{Float64},
    medium::Model, WFSize::WF2DSize; ext = 10, iflag = 2)
    if sn == 1
      if typeof(position) == Array{T1,2} || typeof(origintime) == Array{T2,1} || typeof(sourcetype) == Array{String,1} || typeof(waveform) == Array{Float64,2}
        error("Input source number does not match other parameters. Check if they are multiple sources.")
      else Sour = source(position, origintime, ext, sourcetype, medium, waveform, WFSize, iflag)
      end
    elseif sn > 1
      if typeof(position) == Array{T1,1} || typeof(origintime) == T2 || typeof(sourcetype) == String || typeof(waveform) == Array(Float64,1)
        error("Input source number does not match other parameters. Check if it is a single source.")
      else tmp = MakeMultiSource(sn, sourcetype, position, origintime, waveform)
        Sour = source(sn,tmp,ext,medium,WFSize,iflag)
      end
    end
end

function InitSou{T1<:Real}(sn::Int64, position::Array{T1},
  sourcetype::Union{String,Array{String,1}}, waveform::Array{Float64},
    medium::Model, WFSize::WF2DSize; ext = 10, iflag = 2)
    if sn == 1
      origintime = 0
      if typeof(position) == Array{T1,2} || typeof(sourcetype) == Array{String,1} || typeof(waveform) == Array{Float64,2}
        error("Input source number does not match other parameters. Check if they are multiple sources.")
      else Sour = source(position, origintime, ext, sourcetype, medium, waveform, WFSize, iflag)
      end
    elseif sn > 1
      origintime = zeros(sn)
      if typeof(position) == Array{T1,1} || typeof(sourcetype) == String || typeof(waveform) == Array(Float64,1)
        error("Input source number does not match other parameters. Check if it is a single source.")
      else tmp = MakeMultiSource(sn, sourcetype, position, origintime, waveform)
        Sour = source(sn,tmp,ext,medium,WFSize,iflag)
      end
    end
end

function InitSou{T1,T3<:Real}(sn::Int64, position::Array{T1},
  sourcetype::Union{String,Array{String,1}}, waveform::Array{Float64}, T::T3,
    medium::Model, WFSize::WF2DSize; ext = 10, iflag = 2)
    Tn = Int64(round(T/medium.dt))
    if sn == 1
      origintime = 0
      if typeof(position) == Array{T1,2} || typeof(sourcetype) == Array{String,1} || typeof(waveform) == Array{Float64,2}
        error("Input source number does not match other parameters. Check if they are multiple sources.")
      else Sour = source(position, origintime, ext, sourcetype, medium, waveform,Tn, WFSize, iflag)
      end
    elseif sn > 1
      origintime = zeros(sn)
      if typeof(position) == Array{T1,1} || typeof(sourcetype) == String || typeof(waveform) == Array(Float64,1)
        error("Input source number does not match other parameters. Check if it is a single source.")
      else tmp = MakeMultiSource(sn, sourcetype, position, origintime, waveform)
        Sour = source(sn,tmp,ext,Tn,medium,WFSize,iflag)
      end
    end
end

function InitSou{T1,T2,T3<:Real}(sn::Int64, position::Array{T1}, origintime::Union{T2,Array{T2,1}},
  sourcetype::Union{String,Array{String,1}}, waveform::Array{Float64}, T::T3,
    medium::Model, WFSize::WF2DSize; ext = 10, iflag = 2)
    Tn = Int64(round(T/medium.dt))
    if sn == 1
      if typeof(position) == Array{T1,2} || typeof(origintime) == Array{T2,1} || typeof(sourcetype) == Array{String,1} || typeof(waveform) == Array{Float64,2}
        error("Input source number does not match other parameters. Check if they are multiple sources.")
      else Sour = source(position, origintime, ext, sourcetype, medium, waveform, Tn, WFSize, iflag)
      end
    elseif sn > 1
      if typeof(position) == Array{T1,1} || typeof(origintime) == T2 || typeof(sourcetype) == String || typeof(waveform) == Array(Float64,1)
        error("Input source number does not match other parameters. Check if it is a single source.")
      else tmp = MakeMultiSource(sn, sourcetype, position, origintime, waveform)
        Sour = source(sn,tmp,ext,Tn,medium,WFSize,iflag)
      end
    end
end

function InitSou(rsou::recsou,medium::Model, WFSize::WF2DSize; ext = 10, iflag = 2)

    Tn = rsou.sTn
    position = rsou.slc
    sourcetype = rsou.stp
    waveform = rsou.swf
    sn = rsou.sn
    if sn == 1
      origintime = 0
      Sour = source(position, origintime, ext, sourcetype, medium, waveform, WFSize, iflag)
    elseif sn > 1
      origintime = zeros(sn)
      tmp = MakeMultiSource(sn, sourcetype, position, origintime, waveform)
      Sour = source(sn,tmp,ext,Tn,medium,WFSize,iflag)
    end
end
