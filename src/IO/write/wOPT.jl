#=== write Recording Components ===#
function wCT(fid::IOStream,OptCpnt::Array{String},nwf::nwf2d)

    SVx = nwf.BDnvx[1]*nwf.BDnvx[2];
    SVz = nwf.BDnvz[1]*nwf.BDnvz[2];
    STxx = nwf.BDntpp[1]*nwf.BDntpp[2];
    STzz = nwf.BDntpp[1]*nwf.BDntpp[2];
    STxz = nwf.BDntxz[1]*nwf.BDntxz[2];
    #SWy = nwf.N_BDWy[1]*nwf.N_BDWy[2];

    SizeWFLib = Dict(
    "V"=>SVx+SVz,
    "Vx"=>SVx,
    "Vz"=>SVz,
    "Txx"=>STxx,
    "Tzz"=>STzz,
    "Tii"=>STxx+STzz,
    "Txz"=>STxz,
    #"Wy"=>SWy,
    "full"=>SVx+SVz+STxx+STzz+STxz)#+SWy)

    if "full" in OptCpnt
        Opt = Dict("full"=>"[Vx;Vz;Txx;Tzz;Txz]\n")
        write(fid,Opt["full"])
        NumCpnt = 6;
        SizeWFOpt = SizeWFLib["full"];
    else
        NumCpnt = 0
        SizeWFOpt = 0
        Opt = Dict(
        "V"=>["Vx","Vz"],
        "Tii"=>["Txx","Tzz"])

        for Cpnt in ["V","Tii"]
            if Cpnt in OptCpnt
                filter!(x->!(x in Opt[Cpnt]),OptCpnt)
            end
        end
        
        wOpt = Dict(
        "V"=>"Vx;Vz",
        "Tii"=>"Txx;Tzz")

        write(fid,"[")
        NCpnt = length(OptCpnt)-1==0?1:length(OptCpnt)-1
        for nCpnt in 1:NCpnt
            if OptCpnt[nCpnt] in ["V","Tii"]
                write(fid, wOpt[OptCpnt[nCpnt]])
                NumCpnt = NumCpnt + 2
            else write(fid, OptCpnt[nCpnt])
                NumCpnt = NumCpnt + 1
            end
            if length(OptCpnt) > 1
                write(fid,";")
            end
            SizeWFOpt = SizeWFOpt + SizeWFLib[OptCpnt[nCpnt]]
        end
        if length(OptCpnt) > 1
            if OptCpnt[length(OptCpnt)] in ["V","Tii"]
                write(fid, wOpt[OptCpnt[length(OptCpnt)]])
                NumCpnt = NumCpnt + 2
            else write(fid, OptCpnt[length(OptCpnt)])
                NumCpnt = NumCpnt + 1
            end
            SizeWFOpt = SizeWFOpt + SizeWFLib[OptCpnt[length(OptCpnt)]]
        end
        write(fid,"]\n")
    end
    return NumCpnt, SizeWFOpt
end





#=== write Recording Parameters ===#
function wRPT(fid::IOStream,nwf::nwf2d,BDnrec::Array{Int64,2},
    OptCpnt::Array{String},medium::medium2d)
    write(fid,"recinfo.bin: ") #   Recordings information file
    write(fid,"Simulated recordings information file, structured as\n")
    write(fid,"[Total Receiver Number; Total Time Sampling; Receiver Location Coordinates [Nz Nx] ]\n\n")

    nz, nx = size(BDnrec)
    write(fid,"rec.bin: ") #   Recordings data file
    write(fid,"recordings data file, structured as ")

    NumCpnt,SizeWFOpt = wCT(fid,OptCpnt,nwf)
    write(fid,"File Size $(round(nx*NumCpnt*medium.nT*8/1024/1024,2))MB\n\n")
end





#=== write WaveField Parameters ===#
function wWFPT(fid::IOStream,nwf::nwf2d,
    Slices::Union{StepRange{Int64,Int64},Array{Int64,1}},
    OptCpnt::Array{String})
    write(fid,"wfinfo.bin: ") # Wavefields.bin
    write(fid,"Simulated wavefields information file, structured as\n")
    write(fid,"[Tpp Model Size(PML) [NZ NX]; Total Time Sampling Number; Corresponding Sampling Time ]\n\n")

    write(fid,"wf.bin: ")
    write(fid,"Simulated wavefields data file, structured as ")
    NumCpnt,SizeWFOpt = wCT(fid,OptCpnt,nwf)
    datasize = round(SizeWFOpt*length(Slices)*8/1024/1024/1024,2)
    write(fid,"File Size $(datasize)GB\n\n")
end





#=== write medium parameters ===#
function wOPT(fid::IOStream)
    wSPT(fid)
    write(fid,"No Recordings\n") # Recordings
    write(fid,"No Wavefields\n\n") # Wavefields
end




#=== write medium and recording parameters ===#
function wOPT(
    fid::IOStream,
    nwf::nwf2d,
    BDnrec::Array{Int64,2},
    medium::medium2d,
    OptCpnt::Array{String})

    wSPT(fid)
    wRPT(fid,nwf,BDnrec,OptCpnt,medium)
    write(fid,"No Wavefields\n\n") # Wavefields
end




#=== write medium and wavefield parameters ===#
function wOPT(
    fid::IOStream,
    nwf::nwf2d,
    Slices::Union{StepRange{Int64,Int64},Array{Int64,1}},
    OptCpnt::Array{String})

    wSPT(fid)
    write(fid,"No Recordings\n") # Recordings
    wWFPT(fid,nwf,Slices,OptCpnt)
end




#=== write medium, recording and wavefield parameters ===#
function wOPT(
    fid::IOStream,
    nwf::nwf2d,
    Slices::Union{StepRange{Int64,Int64},Array{Int64,1}},
    BDnrec::Array{Int64,2},
    medium::medium2d,
    OptCpnt::Array{String})

    wSPT(fid)
    wRPT(fid,nwf,BDnrec,OptCpnt,medium)
    wWFPT(fid,nwf,Slices,OptCpnt)
end





#=== write MOR parameters ===#
function wOPT(
    fid::IOStream,
    model::spmod2d)

    wSPT(fid)
    write(fid,"Model order reduction\n") # Recordings

    write(fid,"wfinfo.bin: ") # Wavefields.bin
    write(fid,"Simulated wavefields information file, structured as\n")
    write(fid,"[Tpp Model Size(PML) [NZ;NX];
    Total Time Sampling Number; Corresponding Sampling Time ]\n\n")

    write(fid,"wf.bin: ")
    write(fid,"Simulated wavefields data file, structured as ")
    write(fid,"[Vxx,Vxz,Vzx,Vzz,Txxx,Txxz,Tzzx,Tzzz,Txzx,Txzz]\n")
    NumCpnt,SizeWFOpt = wCT(fid, ["V","Tii","Txz"],model.nwf)
    datasize = round(2*SizeWFOpt*length(Slices)*8/1024/1024/1024,2)
    write(fid,"File Size $(datasize)GB\n\n")

end
