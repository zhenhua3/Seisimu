#=== write MSSB ===#
function wMSSB(WFSize::WF2DSize,Sour::Source,PF::T,
    iflag::Int64, Medium::Model, Path::Union{String, Void},
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
          write(fid,Float64(Medium.dy))
          write(fid,Float64(Medium.dz))
          write(fid,Float64(Medium.dt))
          write(fid,Medium.Tn)
          write(fid,1) # Source number
          write(fid,Float64(PF))
          write(fid,Sour.posn)
          write(fid,Sour.waveform)
          write(fid,Sour.sotn)
          if Sour.stp == "expl"
              write(fid,1)
          elseif Sour.stp == "sfx"
              write(fid,2)
          elseif Sour.stp == "sfy"
              write(fid,3)
          elseif Sour.stp == "sfz"
              write(fid,4)
          elseif Sour.stp == "scx"
              write(fid,5)
          elseif Sour.stp == "scy"
              write(fid,6)
          elseif Sour.stp == "scz"
              write(fid,7)
          elseif Sour.stp == "DC"
              write(fid,8)
          elseif Sour.stp == "wx"
              write(fid,9)
          elseif Sour.stp == "wy"
              write(fid,10)
          elseif Sour.stp == "wz"
              write(fid,11)
          end
          close(fid)
      end
    end
end
