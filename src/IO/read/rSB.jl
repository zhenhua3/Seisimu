function readsou(ParaPath::String)

    fid = open(ParaPath,"r")
    iflag, ext = read(fid,Int64,2)
    modelsize = read(fid,Int64,2) # Model size Nz and Nx with PML boundary
    seek(8*3*modelsize[1]*modelsize[2]+8*2)
    dx,dz,dt = read(fid,Float64,3)
    nT = read(fid,Int64,1)
    pkf = read(fid,Float64,1)
    sn = read(fid,Int64,1)
    nloc = read(fid,Int64,sn*2)
    nloc = reshape(nloc,sn,2)

    loc = [nloc[:,1]*dz nloc[:,2]*dx]
    Waveform = read(fid,Float64,nT*sn)
    Waveform = reshape(Waveform,nT,sn)
    not = read(fid,Int64,sn)
    ot = not.*dt
    tpcode = read(fid,Int64,sn)
    tp = Array{String}(sn)
    for i = 1:sn
      if tpcode[i] == 1
        tp[i] = "expl"
      elseif tpcode[i] == 2
        tp[i] = "sfx"
      elseif tpcode[i] == 3
        tp[i] = "sfz"
      elseif tpcode[i] == 4
        tp[i] = "scx"
      elseif tpcode[i] == 5
        tp[i] = "scz"
      elseif tpcode[i] == 6
        tp[i] = "DC"
      end
    end
    sou = Array{source}[sn]
    for i in 1:sn
        sou[i] = source(tp[i],loc[i,:],nloc[i,:],BDnloc[i,:],ot[i],not[i],waveform[:,i])
    end
    close(fid)
    return sou
end
