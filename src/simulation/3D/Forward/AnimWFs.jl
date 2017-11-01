using PyCall
@pyimport matplotlib.pyplot as plt

function AnimWFs(OptPath::String, WF1::WF2D, WF2::WF2D, WFSize::WF2DSize, medium::Model, source::Array{Source}, FD::NSFDMtx,
  pml::dampCoef, ext::Int64, iflag::Int64,
  sn::Int64, Slices::Union{StepRange{Int64,Int64}, Array{Int64,1}, Int64},
  OptCpnt::Union{String,Array{String,1},Array{String,1}};
  wbox=6, hbox=6, clip = 1.0, cmap="seismic", aspect="auto", interval=40, vmax=1, vmin=-1)
# pathout: the directory which save the wavefield moive

  fig = plt.figure(figsize=(wbox, hbox))
  ims = PyCall.PyObject[]

  if iflag == 1 # free surface
    if OptCpnt == "Vx"
      for it = 1 : medium.Tn
        Addsource!(WF1, sn, source, it)
        Onestep2e!(WF1, WF2, FD, pml)
        if it in Slices
          Snpt = reshape(WF1.VecBDVx, WFSize.N_BDVx[1], WFSize.N_BDVx[2])
          # vmax = maximum(abs(Snpt)) * clip
          # vmin = -vmax
          im = plt.imshow(Snpt[1:WFSize.N_Vx[1],ext+1:ext+WFSize.N_Vx[2]],cmap=cmap,vmin=vmin,vmax=vmax,extent=[medium.dx, WFSize.N_Vx[2]*medium.dx, WFSize.N_Vx[1]*medium.dz, medium.dz],aspect=aspect,clip=clip)
          plt.title("Vx")
          push!(ims, PyCall.PyObject[im])
        end
      end
      ani = anim.ArtistAnimation(fig, ims, interval=interval, blit=true)
      if OptPath[end] == "/"
        pathout = join([OptPath OptCpnt ".mp4"])
        ani[:save](pathout)
      else pathout = join([OptPath "/" OptCpnt ".mp4"])
        ani[:save](pathout)
      end
    elseif OptCpnt == "Vz"
      for it = 1 : medium.Tn
        Addsource!(WF1, sn, source, it)
        Onestep2e!(WF1, WF2, FD, pml)
        if it in Slices
          Snpt = reshape(WF1.VecBDVz, WFSize.N_BDVz[1], WFSize.N_BDVz[2])
          # vmax = maximum(abs(data)) * clip
          # vmin = -vmax
          im = plt.imshow(Snpt[1:WFSize.N_Vz[1],ext+1:ext+WFSize.N_Vz[2]],cmap=cmap,vmin=vmin,vmax=vmax,extent=[medium.dx, WFSize.N_Vz[2]*medium.dx, WFSize.N_Vz[1]*medium.dz, medium.dz],aspect=aspect,clip=clip)
          plt.title("Vz")
          push!(ims, PyCall.PyObject[im])
        end
      end
      ani = anim.ArtistAnimation(fig, ims, interval=interval, blit=true)
      if OptPath[end] == "/"
        pathout = join([OptPath OptCpnt ".mp4"])
        ani[:save](pathout)
      else pathout = join([OptPath "/" OptCpnt ".mp4"])
        ani[:save](pathout)
      end
    elseif OptCpnt == "P"
      for it = 1 : medium.Tn
        Addsource!(WF1, sn, source, it)
        Onestep2e!(WF1, WF2, FD, pml)
        if it in Slices
          Snpt = reshape((WF1.VecBDTxx+WF1.VecBDTzz)/2, WFSize.N_BDTpp[1], WFSize.N_BDTpp[2])
          # vmax = maximum(abs(Snpt)) * clip
          # vmin = -vmax
          im = plt.imshow(Snpt[1:WFSize.N_Tpp[1],ext+1:ext+WFSize.N_Tpp[2]],cmap=cmap,vmin=vmin,vmax=vmax,extent=[medium.dx, WFSize.N_Tpp[2]*medium.dx, WFSize.N_Tpp[1]*medium.dz, medium.dz],aspect=aspect,clip=clip)
          plt.title("P")
          push!(ims, PyCall.PyObject[im])
        end
      end
      ani = anim.ArtistAnimation(fig, ims, interval=interval, blit=true)
      if OptPath[end] == "/"
        pathout = join([OptPath OptCpnt ".mp4"])
        ani[:save](pathout)
      else pathout = join([OptPath "/" OptCpnt ".mp4"])
        ani[:save](pathout)
      end
    end

  elseif iflag == 2 # unlimited medium
    if OptCpnt == "Vx"
      for it = 1 : medium.Tn
        Addsource!(WF1, sn, source, it)
        Onestep2e!(WF1, WF2, FD, pml)
        if it in Slices
          Snpt = reshape(WF1.VecBDVx, WFSize.N_BDVx[1], WFSize.N_BDVx[2])
          # vmax = maximum(abs(Snpt)) * clip
          # vmin = -vmax
          im = plt.imshow(Snpt[ext+1:ext+WFSize.N_Vx[1],ext+1:ext+WFSize.N_Vx[2]],cmap=cmap,vmin=vmin,vmax=vmax,extent=[medium.dx, WFSize.N_Vx[2]*medium.dx, WFSize.N_Vx[1]*medium.dz, medium.dz],aspect=aspect,clip=clip)
          plt.title("Vx")
          push!(ims, PyCall.PyObject[im])
        end
      end
      ani = anim.ArtistAnimation(fig, ims, interval=interval, blit=true)
      if OptPath[end] == "/"
        pathout = join([OptPath OptCpnt ".mp4"])
        ani[:save](pathout)
      else pathout = join([OptPath "/" OptCpnt ".mp4"])
        ani[:save](pathout)
      end
    elseif OptCpnt == "Vz"
      for it = 1 : medium.Tn
        Addsource!(WF1, sn, source, it)
        Onestep2e!(WF1, WF2, FD, pml)
        if it in Slices
          Snpt = reshape(WF1.VecBDVz, WFSize.N_BDVz[1], WFSize.N_BDVz[2])
          # vmax = maximum(abs(data)) * clip
          # vmin = -vmax
          im = plt.imshow(Snpt[ext+1:ext+WFSize.N_Vz[1],ext+1:ext+WFSize.N_Vz[2]],cmap=cmap,vmin=vmin,vmax=vmax,extent=[medium.dx, WFSize.N_Vz[2]*medium.dx, WFSize.N_Vz[1]*medium.dz, medium.dz],aspect=aspect,clip=clip)
          plt.title("Vz")
          push!(ims, PyCall.PyObject[im])
        end
      end
      ani = anim.ArtistAnimation(fig, ims, interval=interval, blit=true)
      if OptPath[end] == "/"
        pathout = join([OptPath OptCpnt ".mp4"])
        ani[:save](pathout)
      else pathout = join([OptPath "/" OptCpnt ".mp4"])
        ani[:save](pathout)
      end
    elseif OptCpnt == "P"
      for it = 1 : medium.Tn
        Addsource!(WF1, sn, source, it)
        Onestep2e!(WF1, WF2, FD, pml)
        if it in Slices
          Snpt = reshape((WF1.VecBDTxx+WF1.VecBDTzz)/2, WFSize.N_BDTpp[1], WFSize.N_BDTpp[2])
          # vmax = maximum(abs(Snpt)) * clip
          # vmin = -vmax
          im = plt.imshow(Snpt[ext+1:ext+WFSize.N_Tpp[1],ext+1:ext+WFSize.N_Tpp[2]],cmap=cmap,vmin=vmin,vmax=vmax,extent=[medium.dx, WFSize.N_Tpp[2]*medium.dx, WFSize.N_Tpp[1]*medium.dz, medium.dz],aspect=aspect,clip=clip)
          plt.title("P")
          push!(ims, PyCall.PyObject[im])
        end
      end
      ani = anim.ArtistAnimation(fig, ims, interval=interval, blit=true)
      if OptPath[end] == "/"
        pathout = join([OptPath OptCpnt ".mp4"])
        ani[:save](pathout)
      else pathout = join([OptPath "/" OptCpnt ".mp4"])
        ani[:save](pathout)
      end
    end
  end
end


function AnimWFs(OptPath::String, WF1::WF2D, WF2::WF2D, WFSize::WF2DSize, medium::Model, source::Source, FD::NSFDMtx,
  pml::dampCoef, ext::Int64, iflag::Int64, Slices::Union{StepRange{Int64,Int64}, Array{Int64,1}, Int64},
  OptCpnt::Union{String,Array{String,1},Array{String,1}};
  wbox=6, hbox=6, clip = 1.0, cmap="seismic", aspect="auto", interval=40, vmax=1, vmin=-1)
# pathout: the directory which save the wavefield moive

  fig = plt.figure(figsize=(wbox, hbox))
  ims = PyCall.PyObject[]

  if iflag == 1 # free surface
    if OptCpnt == "Vx"
      for it = 1 : medium.Tn
        Addsource!(WF1, source, it)
        Onestep2e!(WF1, WF2, FD, pml)
        if it in Slices
          Snpt = reshape(WF1.VecBDVx, WFSize.N_BDVx[1], WFSize.N_BDVx[2])
          # vmax = maximum(abs(Snpt)) * clip
          # vmin = -vmax
          im = plt.imshow(Snpt[1:WFSize.N_Vx[1],ext+1:ext+WFSize.N_Vx[2]],cmap=cmap,vmin=vmin,vmax=vmax,extent=[medium.dx, WFSize.N_Vx[2]*medium.dx, WFSize.N_Vx[1]*medium.dz, medium.dz],aspect=aspect,clip=clip)
          plt.title("Vx")
          push!(ims, PyCall.PyObject[im])
        end
      end
      ani = anim.ArtistAnimation(fig, ims, interval=interval, blit=true)
      if OptPath[end] == "/"
        pathout = join([OptPath OptCpnt ".mp4"])
        ani[:save](pathout)
      else pathout = join([OptPath "/" OptCpnt ".mp4"])
        ani[:save](pathout)
      end
    elseif OptCpnt == "Vz"
      for it = 1 : medium.Tn
        Addsource!(WF1, source, it)
        Onestep2e!(WF1, WF2, FD, pml)
        if it in Slices
          Snpt = reshape(WF1.VecBDVz, WFSize.N_BDVz[1], WFSize.N_BDVz[2])
          # vmax = maximum(abs(data)) * clip
          # vmin = -vmax
          im = plt.imshow(Snpt[1:WFSize.N_Vz[1],ext+1:ext+WFSize.N_Vz[2]],cmap=cmap,vmin=vmin,vmax=vmax,extent=[medium.dx, WFSize.N_Vz[2]*medium.dx, WFSize.N_Vz[1]*medium.dz, medium.dz],aspect=aspect,clip=clip)
          plt.title("Vz")
          push!(ims, PyCall.PyObject[im])
        end
      end
      ani = anim.ArtistAnimation(fig, ims, interval=interval, blit=true)
      if OptPath[end] == "/"
        pathout = join([OptPath OptCpnt ".mp4"])
        ani[:save](pathout)
      else pathout = join([OptPath "/" OptCpnt ".mp4"])
        ani[:save](pathout)
      end
    elseif OptCpnt == "P"
      for it = 1 : medium.Tn
        Addsource!(WF1, source, it)
        Onestep2e!(WF1, WF2, FD, pml)
        if it in Slices
          Snpt = reshape((WF1.VecBDTxx+WF1.VecBDTzz)/2, WFSize.N_BDTpp[1], WFSize.N_BDTpp[2])
          # vmax = maximum(abs(Snpt)) * clip
          # vmin = -vmax
          im = plt.imshow(Snpt[1:WFSize.N_Tpp[1],ext+1:ext+WFSize.N_Tpp[2]],cmap=cmap,vmin=vmin,vmax=vmax,extent=[medium.dx, WFSize.N_Tpp[2]*medium.dx, WFSize.N_Tpp[1]*medium.dz, medium.dz],aspect=aspect,clip=clip)
          plt.title("P")
          push!(ims, PyCall.PyObject[im])
        end
      end
      ani = anim.ArtistAnimation(fig, ims, interval=interval, blit=true)
      if OptPath[end] == "/"
        pathout = join([OptPath OptCpnt ".mp4"])
        ani[:save](pathout)
      else pathout = join([OptPath "/" OptCpnt ".mp4"])
        ani[:save](pathout)
      end
    end

  elseif iflag == 2 # unlimited medium
    if OptCpnt == "Vx"
      for it = 1 : medium.Tn
        Addsource!(WF1, ource, it)
        Onestep2e!(WF1, WF2, FD, pml)
        if it in Slices
          Snpt = reshape(WF1.VecBDVx, WFSize.N_BDVx[1], WFSize.N_BDVx[2])
          # vmax = maximum(abs(Snpt)) * clip
          # vmin = -vmax
          im = plt.imshow(Snpt[ext+1:ext+WFSize.N_Vx[1],ext+1:ext+WFSize.N_Vx[2]],cmap=cmap,vmin=vmin,vmax=vmax,extent=[medium.dx, WFSize.N_Vx[2]*medium.dx, WFSize.N_Vx[1]*medium.dz, medium.dz],aspect=aspect,clip=clip)
          plt.title("Vx")
          push!(ims, PyCall.PyObject[im])
        end
      end
      ani = anim.ArtistAnimation(fig, ims, interval=interval, blit=true)
      if OptPath[end] == "/"
        pathout = join([OptPath OptCpnt ".mp4"])
        ani[:save](pathout)
      else pathout = join([OptPath "/" OptCpnt ".mp4"])
        ani[:save](pathout)
      end
    elseif OptCpnt == "Vz"
      for it = 1 : medium.Tn
        Addsource!(WF1, source, it)
        Onestep2e!(WF1, WF2, FD, pml)
        if it in Slices
          Snpt = reshape(WF1.VecBDVz, WFSize.N_BDVz[1], WFSize.N_BDVz[2])
          # vmax = maximum(abs(data)) * clip
          # vmin = -vmax
          im = plt.imshow(Snpt[ext+1:ext+WFSize.N_Vz[1],ext+1:ext+WFSize.N_Vz[2]],cmap=cmap,vmin=vmin,vmax=vmax,extent=[medium.dx, WFSize.N_Vz[2]*medium.dx, WFSize.N_Vz[1]*medium.dz, medium.dz],aspect=aspect,clip=clip)
          plt.title("Vz")
          push!(ims, PyCall.PyObject[im])
        end
      end
      ani = anim.ArtistAnimation(fig, ims, interval=interval, blit=true)
      if OptPath[end] == "/"
        pathout = join([OptPath OptCpnt ".mp4"])
        ani[:save](pathout)
      else pathout = join([OptPath "/" OptCpnt ".mp4"])
        ani[:save](pathout)
      end
    elseif OptCpnt == "P"
      for it = 1 : medium.Tn
        Addsource!(WF1, source, it)
        Onestep2e!(WF1, WF2, FD, pml)
        if it in Slices
          Snpt = reshape((WF1.VecBDTxx+WF1.VecBDTzz)/2, WFSize.N_BDTpp[1], WFSize.N_BDTpp[2])
          # vmax = maximum(abs(Snpt)) * clip
          # vmin = -vmax
          im = plt.imshow(Snpt[ext+1:ext+WFSize.N_Tpp[1],ext+1:ext+WFSize.N_Tpp[2]],cmap=cmap,vmin=vmin,vmax=vmax,extent=[medium.dx, WFSize.N_Tpp[2]*medium.dx, WFSize.N_Tpp[1]*medium.dz, medium.dz],aspect=aspect,clip=clip)
          plt.title("P")
          push!(ims, PyCall.PyObject[im])
        end
      end
      ani = anim.ArtistAnimation(fig, ims, interval=interval, blit=true)
      if OptPath[end] == "/"
        pathout = join([OptPath OptCpnt ".mp4"])
        ani[:save](pathout)
      else pathout = join([OptPath "/" OptCpnt ".mp4"])
        ani[:save](pathout)
      end
    end
  end
end
