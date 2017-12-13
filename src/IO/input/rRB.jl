type recdata
  component :: Array{String} # provide information of output data component
  recnumber :: Int64 # Number of Receivers
  location :: Array{Float64} # Receiver coordinates
  nlocation :: Array{Int64}
  BDnlocation :: Array{Int64}
  timesamp :: Int64 # Total time sampling Tn
  data :: Array{Float64}
end

function readrec(DataPath::String)
    fid = open(DataPath,"r")
    CpntDict = Dict(1=>"P",2=>"Vx",3=>"Vz",
                    4=>"Txx",5=>"Tzz",6=>"Txz")

        rn = read(fid,Int64,1)
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
        nT = read(fid,Int64,1) # Number of Receivers and Total time sampling Tn
        RowNum = length(component)*rn
        ColNum = nT
        loc = read(fid, Float64, 2*rn)
        loc = reshape(loc, rn,2)
        nloc = read(fid,Int64, 2*rn)
        nloc = reshape(nloc,rn,2)
        BDnloc = read(fid,Int64, 2*rn)
        BDnloc = reshape(BDnloc,rn,2)

        data = read(fid,Float64,RowNum*ColNum)
        data = reshape(data,RowNum,ColNum)

        if eof(fid)
            println("Data reading successful!")
            return recdata(Cpnt,rn,loc,nloc,BDnloc,nT,data)
            close(fid)
        elseif !(eof(fid))
            close(fid)
            error("Check data length, not finished!")
        end

    end
