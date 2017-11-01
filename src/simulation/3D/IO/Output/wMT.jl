#=== write Medium ===#
function wMT(fid::IOStream,WFSize::WF2DSize,Medium::Model)
    write(fid,"Medium\n")
    write(fid,"  Type: ") # free surface or unlimited Medium
    if iflag == 1
      write(fid,"Free surface\n")
    elseif iflag == 2
      write(fid,"Unlimited medium\n")
    end
    write(fid,"  Depth(m): $(Medium.Dep)\n") # Depth
    write(fid,"  HorizonX(m): $(Medium.HorX)\n") # HorizonX
    write(fid,"  HorizonY(m): $(Medium.HorY)\n") # HorizonX
    write(fid,"  Single Side PML Layer Number: $ext\n") # PML
    write(fid,"  Time Step(ms): $(round(Medium.dt*1000,2))\n") # dt
    write(fid,"  Spatial Step(m): $(round(Medium.dx,2))\n") # dx, dy & dz
    write(fid,"  Total Simulation Time(s): $(Medium.TST)\n") # Total Simulation Time
    write(fid,"  Discrete Total Simulation Time: $(Medium.Tn)\n") # Total time sampling number
    write(fid,"  Discrete Model Size(PML) [Z,X,Y]\n") # model size with PML
    write(fid,"     P: [$(WFSize.N_BDTpp[1]),$(WFSize.N_BDTpp[2]),$(WFSize.N_BDTpp[3])]\n")
    write(fid,"   Txz: [$(WFSize.N_BDTxz[1]),$(WFSize.N_BDTxz[2]),$(WFSize.N_BDTxz[3])]\n")
    write(fid,"   Tyz: [$(WFSize.N_BDTyz[1]),$(WFSize.N_BDTyz[2]),$(WFSize.N_BDTyz[3])]\n")
    write(fid,"   Txy: [$(WFSize.N_BDTxz[1]),$(WFSize.N_BDTxz[2]),$(WFSize.N_BDTxy[3])]\n")
    write(fid,"    Vx: [$(WFSize.N_BDVx[1]),$(WFSize.N_BDVx[2]),$(WFSize.N_BDVx[3])]\n")
    write(fid,"    Vy: [$(WFSize.N_BDVy[1]),$(WFSize.N_BDVy[2]),$(WFSize.N_BDVy[3])]\n")
    write(fid,"    Vz: [$(WFSize.N_BDVz[1]),$(WFSize.N_BDVz[2]),$(WFSize.N_BDVz[3])]\n")
    write(fid,"    Wx: [$(WFSize.N_BDWx[1]),$(WFSize.N_BDWx[2]),$(WFSize.N_BDWx[3])]\n")
    write(fid,"    Wy: [$(WFSize.N_BDWy[1]),$(WFSize.N_BDWy[2]),$(WFSize.N_BDWy[3])]\n")
    write(fid,"    Wz: [$(WFSize.N_BDWz[1]),$(WFSize.N_BDWz[2]),$(WFSize.N_BDWz[3])]\n")
    write(fid,"  Discrete Model Size(No PML) [Z,X,Y]\n") # model size with no PML
    write(fid,"     P: [$(WFSize.N_Tpp[1]),$(WFSize.N_Tpp[2]),$(WFSize.N_Tpp[3])]\n")
    write(fid,"   Txz: [$(WFSize.N_Txz[1]),$(WFSize.N_Txz[2]),$(WFSize.N_Txz[3])]\n")
    write(fid,"   Tyz: [$(WFSize.N_Tyz[1]),$(WFSize.N_Tyz[2]),$(WFSize.N_Tyz[3])]\n")
    write(fid,"   Txy: [$(WFSize.N_Txz[1]),$(WFSize.N_Txz[2]),$(WFSize.N_Txy[3])]\n")
    write(fid,"    Vx: [$(WFSize.N_Vx[1]),$(WFSize.N_Vx[2]),$(WFSize.N_Vx[3])]\n")
    write(fid,"    Vy: [$(WFSize.N_Vy[1]),$(WFSize.N_Vy[2]),$(WFSize.N_Vy[3])]\n")
    write(fid,"    Vz: [$(WFSize.N_Vz[1]),$(WFSize.N_Vz[2]),$(WFSize.N_Vz[3])]\n")
    write(fid,"    Wx: [$(WFSize.N_Wx[1]),$(WFSize.N_Wx[2]),$(WFSize.N_Wx[3])]\n")
    write(fid,"    Wy: [$(WFSize.N_Wy[1]),$(WFSize.N_Wy[2]),$(WFSize.N_Wy[3])]\n")
    write(fid,"    Wz: [$(WFSize.N_Wz[1]),$(WFSize.N_Wz[2]),$(WFSize.N_Wz[3])]\n\n")
end
