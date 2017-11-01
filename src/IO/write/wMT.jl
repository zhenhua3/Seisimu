#=== write medium for non split pml ===#
function wMT(fid::IOStream,model::nspmod2d)
    medium = model.medium
    nwf = model.nwf
    iflag = model.medium.iflag

    write(fid,"medium\n")
    write(fid,"  Type: ") # free surface or unlimited medium
    if iflag == 1
      write(fid,"Free surface\n")
    elseif iflag == 2
      write(fid,"Unlimited medium\n")
    end
    write(fid,"  Depth(m): $(medium.DZ)\n") # Depth
    write(fid,"  HorizonX(m): $(medium.HX)\n") # HorizonX
    write(fid,"  Single Side PML Layer Number: $(medium.ext)\n") # PML
    write(fid,"  Time Step(ms): $(round(medium.dt*1000,2))\n") # dt
    write(fid,"  Spatial Step(m): $(round(medium.dx,2))\n") # dx, dy & dz
    write(fid,"  Total Simulation Time(s): $(medium.T)\n") # Total Simulation Time
    write(fid,"  Discrete Total Simulation Time: $(medium.nT)\n") # Total time sampling number
    write(fid,"  Discrete Model Size(PML) [Z,X]\n") # model size with PML
    write(fid,"     P: [$(nwf.BDntpp[1]),$(nwf.BDntpp[2])]\n")
    write(fid,"   Txz: [$(nwf.BDntxz[1]),$(nwf.BDntxz[2])]\n")
    write(fid,"    Vx: [$(nwf.BDnvx[1]),$(nwf.BDnvx[2])]\n")
    write(fid,"    Vz: [$(nwf.BDnvz[1]),$(nwf.BDnvz[2])]\n")
    write(fid,"  Discrete Model Size(No PML) [Z,X,Y]\n") # model size with no PML
    write(fid,"     P: [$(nwf.ntpp[1]),$(nwf.ntpp[2])]\n")
    write(fid,"   Txz: [$(nwf.ntxz[1]),$(nwf.ntxz[2])]\n")
    write(fid,"    Vx: [$(nwf.nvx[1]),$(nwf.nvx[2])]\n")
    write(fid,"    Vz: [$(nwf.nvz[1]),$(nwf.nvz[2])]\n\n")
end

#=== write medium for split pml ===#
function wMT(fid::IOStream,model::spmod2d)
    medium = model.medium
    nwf = model.nwf
    iflag = model.medium.iflag

    write(fid,"medium\n")
    write(fid,"  Type: ") # free surface or unlimited medium
    if iflag == 1
      write(fid,"Free surface\n")
    elseif iflag == 2
      write(fid,"Unlimited medium\n")
    end
    write(fid,"  Depth(m): $(medium.DZ)\n") # Depth
    write(fid,"  HorizonX(m): $(medium.HX)\n") # HorizonX
    write(fid,"  Single Side PML Layer Number: $medium.ext\n") # PML
    write(fid,"  Time Step(ms): $(round(medium.dt*1000,2))\n") # dt
    write(fid,"  Spatial Step(m): $(round(medium.dx,2))\n") # dx, dy & dz
    write(fid,"  Total Simulation Time(s): $(medium.T)\n") # Total Simulation Time
    write(fid,"  Discrete Total Simulation Time: $(medium.nT)\n") # Total time sampling number
    write(fid,"  Discrete Model Size(PML) [Z,X]\n") # model size with PML
    write(fid,"  Tppx: [$(nwf.BDntpp[1]),$(nwf.BDntpp[2])]\n")
    write(fid,"  Tppz: [$(nwf.BDntpp[1]),$(nwf.BDntpp[2])]\n")
    write(fid,"  Txzx: [$(nwf.BDntxz[1]),$(nwf.BDntxz[2])]\n")
    write(fid,"  Txzz: [$(nwf.BDntxz[1]),$(nwf.BDntxz[2])]\n")
    write(fid,"   Vxx: [$(nwf.BDnvx[1]),$(nwf.BDnvx[2])]\n")
    write(fid,"   Vxz: [$(nwf.BDnvx[1]),$(nwf.BDnvx[2])]\n")
    write(fid,"   Vzx: [$(nwf.BDnvz[1]),$(nwf.BDnvz[2])]\n")
    write(fid,"   Vzz: [$(nwf.BDnvz[1]),$(nwf.BDnvz[2])]\n")
    write(fid,"  Discrete Model Size(No PML) [Z,X]\n") # model size with no PML
    write(fid,"  Tppx: [$(nwf.ntpp[1]),$(nwf.ntpp[2])]\n")
    write(fid,"  Tppz: [$(nwf.ntpp[1]),$(nwf.ntpp[2])]\n")
    write(fid,"  Txzx: [$(nwf.ntxz[1]),$(nwf.ntxz[2])]\n")
    write(fid,"  Txzz: [$(nwf.ntxz[1]),$(nwf.ntxz[2])]\n")
    write(fid,"   Vxx: [$(nwf.nvx[1]),$(nwf.nvx[2])]\n")
    write(fid,"   Vxz: [$(nwf.nvx[1]),$(nwf.nvx[2])]\n")
    write(fid,"   Vzx: [$(nwf.nvz[1]),$(nwf.nvz[2])]\n")
    write(fid,"   Vzz: [$(nwf.nvz[1]),$(nwf.nvz[2])]\n\n")
end
