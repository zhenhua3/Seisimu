function simurec(
    model::Union{elmod2d,acmod2d},
    sou::Array{source},
    rec::receiver,
    OptCpnt::Array{String},
    RecDataPath::String)

    nT = model.medium.nT
    fid = open(RecDataPath, "a+")
    OptCpnt = checkrepeat(OptCpnt)

    wIB(fid,rec,model.medium.nT,OptCpnt)
    for it = 1 : nT
        addsou!(model.wf, sou, it)
        extrap2d!(model)
        wDB(fid,model.wf,rec.BDnloc,model.nwf,OptCpnt)
    end
    close(fid)
end
