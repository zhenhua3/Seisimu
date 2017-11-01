#=== write Multiple Sources ===#
function wMST(fid::IOStream,sourcenumber::Int64,Sour::Source,PF::T)
    write(fid,"Source\n") # source
    write(fid,"  Number: $(sourcenumber)\n")
    write(fid,"  Location Coordinates [Z,X,Y] (meter): ")
    for i = 1:sourcenumber
      write(fid," [$(Sour[i].pos[1]),$(Sour[i].pos[2])],$(Sour[i].pos[3])],")
    end
    write(fid,"\n")
    write(fid,"  Type: ")
    for i = 1:sourcenumber
      write(fid," $(Sour[i].stp),")
    end
      write(fid,"\n")
      write(fid,"  Source Origin Time(s): ")
    for i = 1:sourcenumber
      write(fid," $(Sour[i].sot),")
    end
    write(fid,"\n")
    write(fid,"  Discrete Source Origin Time: ")
    for i = 1:sourcenumber
      write(fid," $(Sour[i].sotn),")
    end
    write(fid,"\n")
    write(fid,"  Peak Frequency(Hz): $(PF)\n\n")
end
