#=== write Recording Components ===#
function wCT(OptCpnt::Array{String},WFSize::WF2DSize)
    SVx = WFSize.N_BDVx[1]*WFSize.N_BDVx[2];
    SVz = WFSize.N_BDVz[1]*WFSize.N_BDVz[2];
    SVy = WFSize.N_BDVy[1]*WFSize.N_BDVy[2];
    STxx = WFSize.N_BDTxx[1]*WFSize.N_BDTxx[2];
    STyy = WFSize.N_BDTyy[1]*WFSize.N_BDTyy[2];
    STzz = WFSize.N_BDTzz[1]*WFSize.N_BDTzz[2];
    STxz = WFSize.N_BDTxz[1]*WFSize.N_BDTxz[2];
    STxy = WFSize.N_BDTxy[1]*WFSize.N_BDTxy[2];
    STyz = WFSize.N_BDTyz[1]*WFSize.N_BDTyz[2];
    SWx = WFSize.N_BDWx[1]*WFSize.N_BDWx[2];
    SWy = WFSize.N_BDWy[1]*WFSize.N_BDWy[2];
    SWz = WFSize.N_BDWz[1]*WFSize.N_BDWz[2];
    SizeWFLib = Dict("V"=>SVx+SVy+SVz,"Tii"=>STxx+STyy+STzz,"Tij"=>STxz+STxy+STyz,
    "W"=>SWx+SWy+SWz,"Vx"=>SVx,"Vy"=>SVy,"Vz"=>SVz,"Txx"=>STxx,"Tyy"=>STyy,"Tzz"=>STzz,
    "Txz"=>STxz,"Tyz"=>STyz,"Txy"=>STxy,"Wx"=>SWx,"Wy"=>SWy,"Wz"=>SWz,
    "full"=>SVx+SVy+SVz+STxx+STyy+STzz+STxz+STxy+STyz+SWx+SWy+SWz)
    if "full" in OptCpnt
        Opt = Dict("full"=>"[Vx;Vy;Vz;Txx;Tyy;Tzz;Txz;Tyz;Txz;Wx;Wy;Wz]\n")
        write(fid,Opt["full"])
        NumCpnt = 12;
        SizeWFOpt = SizeWFLib["full"];
    else
        NumCpnt = 0
        SizeWFOpt = 0
        Opt = Dict("V"=>["Vx","Vy","Vz"], "Tii"=>["Txx","Tyy","Tzz"],
        "Tij"=>["Txz","Tyz","Txy"], "W"=>["Wx","Wy","Wz"])
        wOpt = Dict("V"=>"Vx;Vy;Vz", "Tii"=>"Txx;Tyy;Tzz",
        "Tij"=>"Txz;Tyz;Txy", "W"=>"Wx;Wy;Wz")
        Lib = ["V","Tii","Tij","W","Vz","Vx","Vy",
        "Tzz","Txx","Tyy","Txz","Tyz","Txy","Wz","Wx","Wy","P"]
        for Cpnt in ["V","Tii","Tij","W"]
            if Cpnt in OptCpnt
                filter!(x->!(x in Opt[Cpnt]),OptCpnt)
            end
        end
        write(fid,"[")
        NCpnt = length(OptCpnt)-1==0?1:length(OptCpnt)-1
        for nCpnt in 1:NCpnt
            if OptCpnt[nCpnt] in ["V","Tii","Tij","W"]
                write(fid, wOpt[OptCpnt[nCpnt]])
                NumCpnt = NumCpnt + 3
            else write(fid, OptCpnt[nCpnt])
                NumCpnt = NumCpnt + 1
            end
            if NCpnt > 1
                write(fid,";")
            end
            SizeWFOpt = SizeWFOpt + SizeWFLib[OptCpnt[nCpnt]]
        end
        if NCpnt > 1
            if OptCpnt[length(OptCpnt)] in ["V","Tii","Tij","W"]
                write(fid, wOpt[OptCpnt[nCpnt]])
                NumCpnt = NumCpnt + 3
            else write(fid, OptCpnt[nCpnt])
                NumCpnt = NumCpnt + 1
            end
            SizeWFOpt = SizeWFOpt + SizeWFLib[OptCpnt[nCpnt]]
        end
        write(fid,"]\n")
    end
    return NumCpnt, SizeWFOpt
end


#=== write Recording Parameters ===#
function wRPT(fid::IOStream,WFSize::WF2DSize,PML_Rec::Array{Int64,2},
    OptCpnt::Array{String},Medium::Model)
    write(fid,"recinfo.bin: ") #   Recordings information file
    write(fid,"Simulated recordings information file, structured as\n")
    write(fid,"[ Component ID; Total Receiver Number; Total Time Sampling;
    Receiver Location Coordinates [Nz;Nx;Ny] ]\n\n")

    nz, nx = size(PML_Rec)
    write(fid,"rec.bin: ") #   Recordings data file
    write(fid,"recordings data file, structured as ")

    NumCpnt,SizeWFOpt = wCT(OptCpnt,WFSize)
    write(fid,"File Size $(round(nx*NumCpnt*Medium.Tn*8/1024/1024,2))MB\n\n")
end


#=== write WaveField Parameters ===#
function wWFPT(fid::IOStream,WFSize::WF2DSize,
    Slices::Union{StepRange{Int64,Int64},Array{Int64,1}},OptCpnt::Array{String})
    write(fid,"wfinfo.bin: ") # Wavefields.bin
    write(fid,"Simulated wavefields information file, structured as\n")
    write(fid,"[ Component ID; Tpp Model Size(PML) [NZ;NX;NY];
    Total Time Sampling Number; Corresponding Sampling Time ]\n\n")

    write(fid,"wf.bin: ")
    write(fid,"Simulated wavefields data file, structured as ")
    NumCpnt,SizeWFOpt = wCT(OptCpnt,WFSize)
    datasize = round(SizeWFOpt*length(Slices)*8/1024/1024/1024,2)
    write(fid,"File Size $(datasize)GB\n\n")
end


#=== write output parameters ===#
function wOPT(fid::IOStream)
    wSPT(fid)
    write(fid,"No Recordings\n") # Recordings
    write(fid,"No Wavefields\n\n") # Wavefields
end

function wOPT(fid::IOStream,WFSize::WF2DSize,
    PML_Rec::Array{Int64,2},Medium::Model,OptCpnt::Array{String})
    wSPT(fid)
    wRPT(fid,WFSize,PML_Rec,OptCpnt,Medium)
    write(fid,"No Wavefields\n\n") # Wavefields
end

function wOPT(fid::IOStream,WFSize::WF2DSize,
    Slices::Union{StepRange{Int64,Int64},Array{Int64,1}},
    OptCpnt::Array{String})
    wSPT(fid)
    write(fid,"No Recordings\n") # Recordings
    wWFPT(fid,WFSize,Slices,OptCpnt)
end

function wOPT(fid::IOStream,WFSize::WF2DSize,
    Slices::Union{StepRange{Int64,Int64},Array{Int64,1}},
    PML_Rec::Array{Int64,2},Medium::Model,
    OptCpnt::Array{String})
    wSPT(fid)
    wRPT(fid,WFSize,PML_Rec,OptCpnt,Medium)
    wWFPT(fid,WFSize,Slices,OptCpnt)
end
