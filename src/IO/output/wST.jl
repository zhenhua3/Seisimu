#=== write Sources ===#
function wST{T<:Real}(fid::IOStream,sounumber::Int64,
    sou::Array{source},pkf::T)
    write(fid,"Source\n") # source
    write(fid,"  Number: $(sounumber)\n")
    write(fid,"  Location Coordinates [Z,X] (meter): ")
    for i = 1:sounumber
      write(fid," [$(sou[i].loc[1]),$(sou[i].loc[2])],")
    end
    write(fid,"\n")
    write(fid,"  Type: ")
    for i = 1:sounumber
      write(fid," $(sou[i].tp),")
    end
      write(fid,"\n")
      write(fid,"  Source Origin Time(s): ")
    for i = 1:sounumber
      write(fid," $(sou[i].ot),")
    end
    write(fid,"\n")
    write(fid,"  Discrete Source Origin Time: ")
    for i = 1:sounumber
      write(fid," $(sou[i].not),")
    end
    write(fid,"\n")
    write(fid,"  Peak Frequency(Hz): $(pkf)\n\n")
end
