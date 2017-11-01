function simuwf(
    model::nspmod2d,
    sou::Array{source},
    Slices::Union{UnitRange{Int64},StepRange{Int64,Int64},Array{Int64,1},Int64},
    OptCpnt::Array{String},
    WFDataPath::String)

    nT = model.medium.nT
    fid = open(WFDataPath, "a+")

    wIB(fid,model,Slices,OptCpnt)
    for it = 1 : nT
         addsou!(model.wf,sou,it)
         Onestep2D!(model)
         wDB(fid,model.wf,OptCpnt,it,Slices)
    end
    close(fid)
end

#=== model order reduction ===#
function simuwf(
    model::spmod2d,
    source::Array{source},
    Slices::Union{UnitRange{Int64},StepRange{Int64,Int64},Array{Int64,1},Int64},
    WFDataPath::String)

    nT = model.medium.nT
    fid = open(WFDataPath, "a+")

    wIB(fid,model,Slices)
    for it = 1 : nT
         addsou!(model.wf,sou,it)
         Onestep2D!(model)
         wDB(fid,model.wf,it,Slices)
    end
    close(fid)
end
