function simurecwf(
    model::Union{elmod2d,acmod2d},
    source::Array{MTsource},
    rec::receiver,
    RecDataPath::String,
    Slices::Union{UnitRange{Int64},StepRange{Int64,Int64},Array{Int64,1},Int64},
    WFDataPath::String,
    OptCpnt::Array{String},
    )

    nT = model.medium.nT

    #============  write both Recordings and Wavefields ===========#
    fid_wf = open(RecDataPath, "a+")
    fid_record = open(WFDataPath, "a+")
    # wIB(fid_record,rec,model.medium.nT,OptCpnt)
    # wIB(fid_wf,model,Slices,OptCpnt)
    for it = 1 : nT
        for isn in 1:length(sou)
            addmt!(model.wf, sou[isn], it)
        end
          run!(model)
          wDB(fid_record,fid_wf,model.wf,rec.BDnloc,
            model.nwf,OptCpnt,it,Slices)
    end
    close(fid_wf)
end
