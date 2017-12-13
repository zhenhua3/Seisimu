function simuwf(
    model::Union{elmod2d,acmod2d},
    sou::Array{source},
    Slices::Union{UnitRange{Int64},StepRange{Int64,Int64},Array{Int64,1},Int64},
    OptCpnt::Array{String},
    WFDataPath::String)

    nT = model.medium.nT
    fid = open(WFDataPath, "a+")

    wIB(fid,model,Slices,OptCpnt)
    for it = 1 : nT
         addsou!(model.wf,sou,it)
         extrap2d!(model)
         wDB(fid,model.wf,OptCpnt,it,Slices)
    end
    close(fid)
end
