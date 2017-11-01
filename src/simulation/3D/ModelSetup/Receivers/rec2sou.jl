type recsou{T1<:Real}
  sn::Int64
  slc::Array{T1}
  stp::Union{String, Array{String}}
  swf::Array{Float64}
  sTn::Int64
end

function rec2sou(recinfopath::String, recpath::String)

  info = ReadInfo(recinfopath,"rec")
  swf = ReadData(recpath,Float64,info.RowNum,info.ColNum)
  swf = swf'
  Nz, sTn = size(swf)

  if info.Comp == "[Vx; Vz]"
    compnum = 2
    sn = info.RecNum
    sounum = sn*compnum
    stp = Array{String}(sounum)
    stp[1:sn] = "sfx"
    stp[sn+1:2*sn] = "sfz"
    slc = Array{Float64}(2,sounum)
    slc[:,1:sn] = info.RecLoc
    slc[:,sn+1:2*sn] = info.RecLoc
    rsou = recsou(sn,slc,stp,swf,sTn)
  elseif info.Comp == "[P]"
    compnum = 1
    sn = info.RecNum
    sounum = sn*compnum
    if sn == 1
      stp = "expl"
    elseif sn > 1
      stp = Array{String}(sounum)
      stp[1:sn] = "expl"
    end
    slc = info.RecLoc
    rsou = recsou(sn,slc,stp,swf,sTn)
  elseif info.Comp == "[Vx; Vz; P]"
    compnum = 3
    sn = info.RecNum
    sounum = sn*compnum
    stp = Array{String}(sn*compnum)
    stp[1:sn] = "sfx"
    stp[sn+1:2*sn] = "sfz"
    stp[2*sn+1:3*sn] = "expl"
    slc = Array{Float64}(2,sounum)
    slc[:,1:sn] = info.RecLoc
    slc[:,sn+1:2*sn] = info.RecLoc
    slc[:,2*sn+1:3*sn] = info.RecLoc
    rsou = recsou(sounum,slc,stp,swf,sTn)
  end
end

  function rec2trsou(recinfopath::String, recpath::String)

    info = ReadInfo(recinfopath,"rec")


    if info.Comp == "[Vx; Vz]"
      compnum = 2
      sn = info.RecNum
      swf = ReadData(recpath,Float64,info.RowNum,info.ColNum)
      swf[:,1:end] = swf[:,end:-1:1]
      swf = swf'
      sTn, Nz = size(swf)
      sounum = sn*compnum
      stp = Array{String}(sounum)
      stp[1:sn] = "sfx"
      stp[sn+1:2*sn] = "sfz"
      slc = Array{Float64}(2,sounum)
      slc[:,1:sn] = info.RecLoc
      slc[:,sn+1:2*sn] = info.RecLoc
      rsou = recsou(sn,slc,stp,swf,sTn)
    elseif info.Comp == "[P]"
      compnum = 1
      sn = info.RecNum
      swf = ReadData(recpath,Float64,info.RowNum,info.ColNum)
      if sn == 1
        swf[1:end] = swf[end:-1:1]
      else
        swf[:,1:end] = swf[:,end:-1:1]
      end
      swf = swf'
      sTn, Nz = size(swf)
      sounum = sn*compnum
      if sn == 1
        stp = "expl"
      elseif sn > 1
        stp = Array{String}(sounum)
        stp[1:sn] = "expl"
      end
      slc = info.RecLoc
      rsou = recsou(sn,slc,stp,swf,sTn)
    elseif info.Comp == "[Vx; Vz; P]"
      compnum = 3
      sn = info.RecNum
      swf = ReadData(recpath,Float64,info.RowNum,info.ColNum)
      swf[:,1:end] = swf[:,end:-1:1]
      swf = swf'
      sTn, Nz = size(swf)
      sounum = sn*compnum
      stp = Array{String}(sn*compnum)
      stp[1:sn] = "sfx"
      stp[sn+1:2*sn] = "sfz"
      stp[2*sn+1:3*sn] = "expl"
      slc = Array{Float64}(2,sounum)
      slc[:,1:sn] = info.RecLoc
      slc[:,sn+1:2*sn] = info.RecLoc
      slc[:,2*sn+1:3*sn] = info.RecLoc
      rsou = recsou(sounum,slc,stp,swf,sTn)
    end
end
