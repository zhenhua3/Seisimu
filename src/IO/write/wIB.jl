#=== rec ===#
function wIB(
    fid::IOStream,
    rec::receiver,
    nT::Int64,
    OptCpnt::Array{String})

    CpntDict = Dict("P"=>1,"Vx"=>2,"Vz"=>3,
                "Txx"=>4,"Tzz"=>5,"Txz"=>6,"V"=>23,"Tii"=>45)


    write(fid,rec.rn) # Total Receiver Number
    write(fid,length(OptCpnt))
    for cpnt in OptCpnt
        write(fid,CpntDict[cpnt])
    end
    write(fid,nT) # Total Time Sampling Tn
    write(fid,rec.loc) #physical location
    write(fid,rec.nloc) #discrete location without PML
    write(fid,rec.BDnloc) #discrete location with PML
end

#=== wavefield ===#
function wIB(
    fid::IOStream,
    model::Union{nspmod2d},
    Slices::Union{UnitRange{Int64},StepRange{Int64,Int64},Array{Int64,1},Int64},
    OptCpnt::Array{String})

    CpntDict = Dict("P"=>1,"Vx"=>2,"Vz"=>3,
                "Txx"=>4,"Tzz"=>5,"Txz"=>6,
                "V"=>23,"Tii"=>45)

    write(fid,model.medium.iflag)
    write(fid,length(OptCpnt))
    for cpnt in OptCpnt
        write(fid,CpntDict[cpnt])
    end
    write(fid,model.nwf.BDntpp)
    write(fid,model.medium.ext)
    write(fid,model.medium.dx)
    write(fid,model.medium.dz)
    write(fid,length(Slices))
    write(fid,Slices.*model.medium.dt)

end

#== mor ===#
function wIB(
    fid::IOStream,
    model::Union{spmod2d},
    Slices::Union{UnitRange{Int64},StepRange{Int64,Int64},Array{Int64,1},Int64})

    write(fid,model.medium.iflag)
    write(fid,model.nwf.BDntpp)
    write(fid,model.medium.dx)
    write(fid,model.medium.dz)
    write(fid,length(Slices))
    write(fid,Slices.*model.medium.dt)

end
