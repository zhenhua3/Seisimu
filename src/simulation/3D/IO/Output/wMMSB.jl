#=== write MMSB ===#
function wMMSB(WFSize::WF2DSize,Sour::Array{Source},sourcenumber::Int64,
    PF::T,iflag::Int64,Medium::Model,Path::Union{String, Void},
    Path_Medium::Union{String, Void},ext::Int64)
    if Path_Medium != nothing
        if Path_Medium[end-2:end] != "bin"
            error("Check output file extensions. It has to be in 'bin' format.")
        else
            fid = open(Path_Medium,"a+")
            write(fid,iflag)
            write(fid,ext)
            write(fid,WFSize.N_BDTpp)
            write(fid,Medium.VP)
            write(fid,Medium.VS)
            write(fid,Medium.Rho)
            write(fid,Float64(Medium.dx))
            write(fid,Float64(Medium.dz))
            write(fid,Float64(Medium.dt))
            write(fid,Medium.Tn)
            write(fid,sourcenumber)
            write(fid,Float64(PF))
            for i = 1:sourcenumber
                write(fid,Sour[i].posn)
            end
            for i = 1:sourcenumber
                write(fid,Sour[i].waveform)
            end
            for i = 1:sourcenumber
                write(fid,Sour[i].sotn)
            end
            for i = 1:sourcenumber
                if Sour[i].stp == "expl"
                    write(fid,1)
                elseif Sour[i].stp == "sfx"
                    write(fid,2)
                elseif Sour[i].stp == "sfy"
                    write(fid,3)
                elseif Sour[i].stp == "sfz"
                    write(fid,4)
                elseif Sour[i].stp == "scx"
                    write(fid,5)
                elseif Sour[i].stp == "scy"
                    write(fid,6)
                elseif Sour[i].stp == "scz"
                    write(fid,7)
                elseif Sour[i].stp == "DC"
                    write(fid,8)
                elseif Sour[i].stp == "wx"
                    write(fid,9)
                elseif Sour[i].stp == "wy"
                    write(fid,10)
                elseif Sour[i].stp == "wz"
                    write(fid,11)
                end
            end
            close(fid)
        end
    end
end
