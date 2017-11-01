type wfinfo
    mediumtype::String
    component :: Array{String} # provide information of output data component
    mediumsize :: Array{Int64}
    ext :: Int64
    gridsize :: Array{Float64}
    timesamp :: Int64 # Total time sampling nT
    timepoint :: Array{Float64}
end

function readwfinfo(DataPath::String)
    fid = open(DataPath,"r")
    iflag = read(fid,Int64,1)
    if iflag == 1
        mediumtype = "free surface"
    elseif iflag == 2
        mediumtype = "unlimited"
    end

    CpntDict = Dict(1=>"P",2=>"Vx",3=>"Vz",
                  4=>"Txx",5=>"Tzz",6=>"Txz")
    nCpnt = read(fid,Int64,1)
    NumCpnt = read(fid,Int64,nCpnt)
    component = []
    for i in NumCpnt
      if i == 23
          push!(component,CpntDict[2])
          push!(component,CpntDict[3])
      elseif i == 45
          push!(component,CpntDict[4])
          push!(component,CpntDict[5])
      else push!(component,CpntDict[i])
      end
    end
    mediumsize = read(fid,Int64,2)
    ext = read(fid,Int64,1)
    gridsize = read(fid,Float64,2)
    timesamp = read(fid,Int64,1)
    timepoint = read(fid,Float64,timesamp)

    return wfinfo(mediumtype, component, mediumsize,
    ext, gridsize, timesamp, timepoint)

end
function readwfdata(fid::IOStream,compnumber::Int64,timesamp::Int64,
    startpoint::Int64,length::Int64)
    # seek(fid,8)
    # nCpnt = read(fid,Int64,1)
    # seek(fid,8*(7+nCpnt))
    # timesamp = seek(fid,Int64,1)
    seek(fid,8*(8+compnumber+timesamp+startpoint))
    requireddata = read(fid,Float64,length)
    return requireddata
end

function readwfdata(fid::IOStream,startpoint::Int64,length::Int64)
    seek(fid,8)
    compnumber = read(fid,Int64,1)
    seek(fid,8*(7+compnumber))
    timesamp = seek(fid,Int64,1)
    seek(fid,8*(8+compnumber+timesamp+startpoint))
    requireddata = read(fid,Float64,length)
    return(requireddata)
end

function readwfdata(fid::IOStream,length::Int64)
    requireddata = read(fid,Float64,length)
    return requireddata 
end

function readwfdata(DataPath::String)
    fid = open(DataPath,"r")
    iflag = read(fid,Int64,1)
    nCpnt = read(fid,Int64,1)
    NumCpnt = read(fid,Int64,nCpnt)
    mediumsize = read(fid,Int64,2)
    ext = read(fid,Int64,1)
    gridsize = read(fid,Float64,2)
    timesamp = read(fid,Int64,1)
    timepoint = read(fid,Float64,timesamp)
    CpntDict = Dict(1=>mediumsize[1]*mediumsize[2],
                    2=>mediumsize[1]*(mediumsize[2]-1),
                    3=>(mediumsize[1]-1)*mediumsize[2],
                    4=>mediumsize[1]*mediumsize[2],
                    5=>mediumsize[1]*mediumsize[2],
                    6=>(mediumsize[1]-1)*(mediumsize[2]-1))
    totalCpnt = 0
    for Cpnt in NumCpnt
        if Cpnt == 23
            totalCpnt = totalCpnt + CpntDict[2] + CpntDict[3]
        elseif Cpnt == 45
            totalCpnt = totalCpnt + CpntDict[4] + CpntDict[5]
        else totalCpnt = totalCpnt + CpntDict[Cpnt]
        end
    end
    requireddata = read(fid,Float64,totalCpnt*timesamp)
    return requireddata
end
