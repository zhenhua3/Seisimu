function checkrepeat(OptCpnt)
    Opt = Dict(
    "V"=>["Vx","Vz"],
    "Tii"=>["Txx","Tzz"])

    for Cpnt in ["V","Tii"]
        if Cpnt in OptCpnt
            filter!(x->!(x in Opt[Cpnt]),OptCpnt)
        end
    end
    return OptCpnt
end

#=== write Recordings into bin file ===#
function wDB(
    fid::IOStream,
    wf::nspwf2d,
    BDnloc::Array{Int64,2},
    nwf::nwf2d,
    OptCpnt::Array{String})

    RecLocLib = Dict(
    "P"=>(BDnloc[2,:]-1).*nwf.BDntpp[1] + BDnloc[1,:],
    "Vx"=>(BDnloc[2,:]-1).*nwf.BDnvx[1] + BDnloc[1,:],
    "Vz"=>(BDnloc[2,:]-1).*nwf.BDnvz[1] + BDnloc[1,:],
    "Txx"=>(BDnloc[2,:]-1).*nwf.BDntpp[1] + BDnloc[1,:],
    "Tzz"=>(BDnloc[2,:]-1).*nwf.BDntpp[1] + BDnloc[1,:],
    "Txz"=>(BDnloc[2,:]-1).*nwf.BDntxz[1] + BDnloc[1,:])
    #"Wy"=>(BDnloc[2,:]-1).*nwf.N_BDWy[1] + BDnloc[1,:])

    RecLib = Dict(
    "P"=>vec((wf.txx+wf.tzz)/2)[RecLocLib["P"]],
    "Vx"=>vec(wf.vx)[RecLocLib["Vx"]],
    "Vz"=>vec(wf.vz)[RecLocLib["Vz"]],
    "Txx"=>vec(wf.txx)[RecLocLib["Txx"]],
    "Txz"=>vec(wf.txz)[RecLocLib["Txz"]],
    "Tzz"=>vec(wf.tzz)[RecLocLib["Tzz"]],
    "V"=>[vec(wf.vx)[RecLocLib["Vx"]];
    vec(wf.vz)[RecLocLib["Vz"]]],
    "Tii"=>[vec(wf.txx)[RecLocLib["Txx"]];
    vec(wf.tzz)[RecLocLib["Tzz"]]])
    #"Wy"=>wf.VecBDWy)

    for OptC in OptCpnt
        write(fid,RecLib[OptC])
    end
end


#=== write Wavefields into bin file ===#
function wDB(
    fid::IOStream,
    wf::nspwf2d,
    OptCpnt::Array{String},
    it::Int64,
    Slices::Union{UnitRange{Int64},StepRange{Int64,Int64},Array{Int64,1},Int64})

    WFLib = Dict(
    "P"=>vec((wf.txx+wf.tzz)/2),
    "Vx"=>vec(wf.vx),
    "Vz"=>vec(wf.vz),
    "Txx"=>vec(wf.txx),
    "Tzz"=>vec(wf.tzz),
    "Txz"=>vec(wf.txz),
    "V"=>[vec(wf.vx);vec(wf.vz)],
    "Tii"=>[vec(wf.txx);vec(wf.tzz)])
    #"Wy"=>wf.VecBDWy)
    if it in Slices
        for OptC in OptCpnt
            write(fid,DataLib[OptC])
        end
    end
end

#=== write both Recordings and Wavefields into bin file ===#
function wDB(
    fidRec::IOStream,
    fidWF::IOStream,
    wf::nspwf2d,
    BDnloc::Array{Int64,2},
    nwf::nwf2d,
    OptCpnt::Array{String},
    it::Int64,
    Slices::Union{UnitRange{Int64},StepRange{Int64,Int64},Array{Int64,1},Int64})

    RecLocLib = Dict(
    "P"=>(BDnloc[2,:]-1).*nwf.BDntpp[1] + BDnloc[1,:],
    "Vx"=>(BDnloc[2,:]-1).*nwf.BDnvx[1] + BDnloc[1,:],
    "Vz"=>(BDnloc[2,:]-1).*nwf.BDnvz[1] + BDnloc[1,:],
    "Txx"=>(BDnloc[2,:]-1).*nwf.BDntpp[1] + BDnloc[1,:],
    "Tzz"=>(BDnloc[2,:]-1).*nwf.BDntpp[1] + BDnloc[1,:],
    "Txz"=>(BDnloc[2,:]-1).*nwf.BDntxz[1] + BDnloc[1,:])
    #"Wy"=>(BDnloc[2,:]-1).*nwf.N_BDWy[1] + BDnloc[1,:])

    RecLib = Dict(
    "P"=>vec((wf.txx+wf.tzz)/2)[RecLocLib["P"]],
    "Vx"=>vec(wf.vx)[RecLocLib["Vx"]],
    "Vz"=>vec(wf.vz)[RecLocLib["Vz"]],
    "Txx"=>vec(wf.txx)[RecLocLib["Txx"]],
    "Txz"=>vec(wf.txz)[RecLocLib["Txz"]],
    "Tzz"=>vec(wf.tzz)[RecLocLib["Tzz"]],
    "V"=>[vec(wf.vx)[RecLocLib["Vx"]];
    vec(wf.vz)[RecLocLib["Vz"]]],
    "Tii"=>[vec(wf.txx)[RecLocLib["Txx"]];
    vec(wf.tzz)[RecLocLib["Tzz"]]])
    #"Wy"=>wf.VecBDWy)

    WFLib = Dict(
    "P"=>vec((wf.txx+wf.tzz)/2),
    "Vx"=>vec(wf.vx),
    "Vz"=>vec(wf.vz),
    "Txx"=>vec(wf.txx),
    "Tzz"=>vec(wf.tzz),
    "Txz"=>vec(wf.txz),
    "V"=>[vec(vx);vec(wf.vz)],
    "Tii"=>[vec(txx);vec(tzz)])
    #"Wy"=>wf.VecBDWy)

    for OptC in OptCpnt
        write(fidRec,RecLib[OptC])
        if it in Slices
            write(fid,WFLib[OptC])
        end
    end

end

#=== write Wavefields into bin file ===#
function wDB(
    fid::IOStream,
    wf::spwf2d,
    it::Int64,
    Slices::Union{UnitRange{Int64},StepRange{Int64,Int64},Array{Int64,1},Int64})

    DataLib = [ vec(wf.vxx);  vec(wf.vxz);  vec(wf.vzx);  vec(wf.vzz);
               vec(wf.txxx); vec(wf.txxz); vec(wf.tzzx); vec(wf.tzzz);
               vec(wf.txzx); vec(wf.txzz)]

    if it in Slices
        write(fid,DataLib)
    end
end
