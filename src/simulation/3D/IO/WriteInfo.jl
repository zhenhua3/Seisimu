#=== No.1 Single/Multiple Source ===#
#=== Output : MSST, MSSB /  MMST, MMSB with no RC and WF ===#
function Writeinfo{T<:Real}(WFSize::WF2DSize,Sour::Source,sourcenumber::Int64,
    PF::T,iflag::Int64, Medium::Model, Path::Union{String, Void},
    Path_Medium::Union{String, Void},ext::Int64)

#=== Simu&OutputInfo.txt ===#
    if Path != nothing
        if Path[end-2:end] != "txt"
            error("Check file extensions. It has to be in 'txt' format. ")
        else
            fid = open(Path,"a+")
            wMT(fid,WFSize,Medium)
            sourcenumber==1?wSST(fid,Sour,PF):wMST(fid,sourcenumber,Sour,PF)
            wSPT(fid)
            wOPT(fid)
            close(fid)
        end
    end
#=== SimuPara.bin ===#
    sourcenumber==1?wMSSB(WFSize, Sour, PF, iflag, Medium, Path, Path_Medium, ext):
    wMMSB(WFSize, Sour, sourcenumber, PF, iflag, Medium, Path, Path_Medium, ext)
end
#======#

#=== No.2 Single/Multiple Source ===#
#=== Output: MSST, MSSB / MMST, MMSB with RC and no WF ===#
function Writeinfo{T<:Real}(WFSize::WF2DSize,Sour::Source,sourcenumber::Int64,
    PF::T,iflag::Int64,Medium::Model,PML_Rec::Array{Int64,2},Path::Union{String, Void},
    Path_Medium::Union{String, Void}, Output::Array{String},ext::Int64)

#=================== Simu&OutputInfo.txt ====================#
    if Path != nothing
        if Path[end-2:end] != "txt"
            error("Check file extensions. It has to be in 'txt' format. ")
        else
            fid = open(Path,"a+")

            wMT(fid,WFSize,Medium)
            sourcenumber==1?wSST(fid,Sour,PF):wMST(fid,sourcenumber,Sour,PF)
            wSPT(fid)
            wOPT(fid,PML_Rec,Medium)
            close(fid)
        end
    end
#=== SimuPara.bin ===#
sourcenumber==1?wMSSB(WFSize, Sour, PF, iflag, Medium, Path, Path_Medium, ext):
wMMSB(WFSize, Sour, sourcenumber, PF, iflag, Medium, Path, Path_Medium, ext)
end
#======#

#=== No.3 Single/Multiple source===#
#=== Output : MSST, MSSB / MMST, MMSB with WF and no RC ===#
function Writeinfo{T<:Real}(WFSize::WF2DSize,Sour::Source,sourcenumber::Int64,
    PF::T,iflag::Int64,Medium::Model,Slices::Union{StepRange{Int64,Int64},Array{Int64,1}},
    Path::Union{String, Void},Path_Medium::Union{String, Void}, Output::Array{String},ext::Int64)

#=================== Simu&OutputInfo.txt ====================#
    if Path != nothing
        if Path[end-2:end] != "txt"
            error("Check file extensions. It has to be in 'txt' format. ")
        else
            fid = open(Path,"a+")

            wMT(fid,WFSize,Medium)
            sourcenumber==1?wSST(fid,Sour,PF):wMST(fid,sourcenumber,Sour,PF)
            wSPT(fid)
            wOPT(fid,WFSize,Slices,OptCpnt)
            close(fid)
        end
    end
#=================== SimuPara.bin ====================#
    sourcenumber==1?wMSSB(WFSize, Sour, PF, iflag, Medium, Path, Path_Medium, ext):
    wMMSB(WFSize, Sour, sourcenumber, PF, iflag, Medium, Path, Path_Medium, ext)
end
#======#

#=== No.4 Single/Multiple source===#
#=== Output : MSST, MSSB / MMST, MMSB with WF and RC ===#
function Writeinfo{T<:Real}(WFSize::WF2DSize,
    Sour::Source,sourcenumber::Int64,PF::T,iflag::Int64,
    Medium::Model,PML_Rec::Array{Int64,2},
    Slices::Union{StepRange{Int64,Int64},Array{Int64,1}},
    Path::Union{String, Void}, Path_Medium::Union{String, Void},
    Output::Array{String},ext::Int64)

#=================== Simu&OutputInfo.txt ====================#
    if Path != nothing
        if Path[end-2:end] != "txt"
            error("Check file extensions. It has to be in 'txt' format. ")
        else
            fid = open(Path,"a+")
            wMT(fid,WFSize,Medium)
            sourcenumber==1?wSST(fid,Sour,PF):wMST(fid,sourcenumber,Sour,PF)
            wSPT(fid)
            wOPT(fid,WFSize,Slices,PML_Rec,Medium,OptCpnt)
            close(fid)
        end
    end
#=================== SimuPara.bin ====================#
    sourcenumber==1?wMSSB(WFSize, Sour, PF, iflag, Medium, Path, Path_Medium, ext):
    wMMSB(WFSize, Sour, sourcenumber, PF, iflag, Medium, Path, Path_Medium, ext)
end
#======#

#=== No.5 ===#
#=== write Single/Multiple source simulation information for Model Order Reduction ===#
function Writeinfo{T<:Real}(WFSize::WF2DSize,Sour::Source,
    PF::T,iflag::Int64,Medium::Model,
    Slices::Union{StepRange{Int64,Int64},Array{Int64,1}},Path::Union{String, Void},
    Path_Medium::Union{String, Void},ext::Int64)

#=================== Simu&OutputInfo.txt ====================#
    if Path != nothing
        if Path[end-2:end] != "txt"
            error("Check file extensions. It has to be in 'txt' format. ")
        else
            fid = open(Path,"a+")

            wMT(fid,WFSize,Medium)
            sourcenumber==1?wSST(fid,Sour,PF):wMST(fid,sourcenumber,Sour,PF)
            wSPT(fid)

            write(fid,"wfinfo.bin: ") # Wavefields.bin
            write(fid,"Simulated wavefields information file, structured as\n")
            write(fid,"[ Component ID; Tpp Model Size(PML) [NZ;Nx]; Total Time Sampling Number; Corresponding Sampling Time ]\n\n")

            write(fid,"wf.bin: ")
            write(fid,"Simulated wavefields data file, structured as ")
            write(fid,"[Vxx,Vxz,Vzx,Vzz,Txxx,Txxz,Tzzx,Tzzz,Txzx,Txzz]\n")
            datasize = 2*round((WFSize.N_BDVx[1]*WFSize.N_BDVx[2]+WFSize.N_BDVz[1]*WFSize.N_BDVz[2]+WFSize.N_BDTpp[1]*WFSize.N_BDTpp[2]+WFSize.N_BDTpp[1]*WFSize.N_BDTpp[2]+WFSize.N_BDTxz[1]*WFSize.N_BDTxz[2])*length(Slices)*8/1024/1024/1024,2)
            write(fid,"File Size $(datasize)GB\n\n")

            write(fid,"No Recordings\n\n") # Recordings
            close(fid)
        end
    end
#=================== SimuPara.bin ====================#
    sourcenumber==1?wMSSB(WFSize, Sour, PF, iflag, Medium, Path, Path_Medium, ext):
    wMMSB(WFSize, Sour, sourcenumber, PF, iflag, Medium, Path, Path_Medium, ext)
end
#======#
