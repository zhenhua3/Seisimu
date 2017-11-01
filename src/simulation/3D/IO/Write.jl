function WriteBasis{T<:Real}(data::Array{T,2}, pathout::String)

  Nz, Nx = size(data)
  fid = open(pathout,"w+")
  write(fid,Nz)
  write(fid,Nx)
  write(fid,data)
  close(fid)

end
