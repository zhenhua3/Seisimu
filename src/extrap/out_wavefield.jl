function simuwf(
    model::Union{elmod2d,acmod2d},
    sou::Array{MTsource},
    Slices::Union{UnitRange{Int64},StepRange{Int64,Int64},Array{Int64,1},Int64},
    OptCpnt::Array{String},
    WFDataPath::String)

    nT = model.medium.nT
    fid = open(WFDataPath, "a+")

    # wIB(fid,model,Slices,OptCpnt)
    for it = 1 : nT
        for isn in 1:length(sou)
            addmt!(model.wf, sou[isn], it)
        end
         run!(model)
         wDB(fid,model.wf,OptCpnt,it,Slices)
    end
    close(fid)
end
