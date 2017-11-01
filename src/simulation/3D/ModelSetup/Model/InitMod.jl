
############ homogeneous acoustic medium ############
function InitMod{T1,T2,T3,T4,T5,T6<:Real}(v_p::T1,
  rho::T2,Depth::T3,Horizon::T4,PeakFreq::T5,T::T6; ext = 10, iflag = 2, mortraining = false, morforward = false)

  medium = model(v_p, rho, Depth, Horizon, PeakFreq, T, ext, iflag)
  WFSize = InitWF2DSize(medium.NDep, medium.NHor, ext, iflag)
  ############### Finite Difference Matrix ################
  FDC = FDCoeff(4)
# =================== unsplit PML boundary ===================================#
  if mortraining == false && morforward == false
    #%%%%%%%% wavefield initialization %%%%%%%%%
    WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
    WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)

    FD = Fdmtx(WFSize, medium, FDC, ext)
    pml = DampBound(WFSize, medium, ext, iflag)
    return medium, WF1, WF2, WFSize, FD, pml
# =================== split PML boundary ===================================#
  elseif mortraining == true && morforward == false
    WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
    WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)

    FD = Fdmtx(WFSize, medium, FDC, ext, iflag)
    return medium, WF1, WF2, WFSize, FD
#=================== MOR forward ============================#
  elseif morforward == true
    FD = Fdmtx(WFSize, medium, FDC, ext, iflag)
    return medium, WFSize, FD
  end

end

############ homogeneous elastic medium ############
function InitMod{T1,T2,T3,T4,T5,T6,T7<:Real}(v_p::T1,
  v_s::T2, rho::T3, Depth::T4, Horizon::T5, PeakFreq::T6, T::T7; ext = 10, iflag = 2, mortraining = false, morforward = false)

  medium = model(v_p, v_s, rho, Depth, Horizon, PeakFreq, T, ext, iflag)
  #%%%%%%%% wavefield initialization %%%%%%%%%
  WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
  WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
  WFSize = InitWF2DSize(medium.NDep, medium.NHor, ext, iflag)
  ############### Finite Difference Matrix ################
  FDC = FDCoeff(4)

  # =================== unsplit PML boundary ===================================#
    if mortraining == false && morforward == false
      #%%%%%%%% wavefield initialization %%%%%%%%%
      WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
      WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)

      FD = Fdmtx(WFSize, medium, FDC, ext)
      pml = DampBound(WFSize, medium, ext, iflag)
      return medium, WF1, WF2, WFSize, FD, pml
  # =================== split PML boundary ===================================#
    elseif mortraining == true && morforward == false
      WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
      WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)

      FD = Fdmtx(WFSize, medium, FDC, ext, iflag)
      return medium, WF1, WF2, WFSize, FD
  #=================== MOR forward ============================#
    elseif morforward == true
      FD = Fdmtx(WFSize, medium, FDC, ext, iflag)
      return medium, WFSize, FD
    end
end

############ elastic layered medium with constant density ############
function InitMod{T1,T2,T3,T4,T5,T6<:Real}(Max_Min_v::Array{T1,1},
  rho::T2, LayerNo::Int64,Depth::T3, Horizon::T4, PeakFreq::T5, T::T6; ext = 10, iflag = 2, mortraining = false, morforward = false)

  medium = model(Max_Min_v, rho, LayerNo, Depth, Horizon, PeakFreq, T, ext, iflag)
  #%%%%%%%% wavefield initialization %%%%%%%%%
  WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
  WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
  WFSize = InitWF2DSize(medium.NDep, medium.NHor, ext, iflag)
  ############### Finite Difference Matrix ################
  FDC = FDCoeff(4)

  # =================== unsplit PML boundary ===================================#
    if mortraining == false && morforward == false
      #%%%%%%%% wavefield initialization %%%%%%%%%
      WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
      WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)

      FD = Fdmtx(WFSize, medium, FDC, ext)
      pml = DampBound(WFSize, medium, ext, iflag)
      return medium, WF1, WF2, WFSize, FD, pml
  # =================== split PML boundary ===================================#
    elseif mortraining == true && morforward == false
      WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
      WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)

      FD = Fdmtx(WFSize, medium, FDC, ext, iflag)
      return medium, WF1, WF2, WFSize, FD
  #=================== MOR forward ============================#
    elseif morforward == true
      FD = Fdmtx(WFSize, medium, FDC, ext, iflag)
      return medium, WFSize, FD
    end
end

############ elastic layered medium with layerd density and constant layer depth ############
function InitMod{T1,T2,T3,T4,T5,T6<:Real}(Max_Min_v::Array{T1,1},
  Max_Min_rho::Array{T2,1}, LayerNo::Int64, Depth::T3, Horizon::T4, PeakFreq::T5, T::T6; ext = 10, iflag = 2, mortraining = false, morforward = false)

  medium = model(Max_Min_v, Max_Min_rho, LayerNo, Depth, Horizon, PeakFreq, T, ext, iflag)
  #%%%%%%%% wavefield initialization %%%%%%%%%
  WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
  WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
  WFSize = InitWF2DSize(medium.NDep, medium.NHor, ext, iflag)
  ############### Finite Difference Matrix ################
  FDC = FDCoeff(4)

  # =================== unsplit PML boundary ===================================#
    if mortraining == false && morforward == false
      #%%%%%%%% wavefield initialization %%%%%%%%%
      WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
      WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)

      FD = Fdmtx(WFSize, medium, FDC, ext)
      pml = DampBound(WFSize, medium, ext, iflag)
      return medium, WF1, WF2, WFSize, FD, pml
  # =================== split PML boundary ===================================#
    elseif mortraining == true && morforward == false
      WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
      WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)

      FD = Fdmtx(WFSize, medium, FDC, ext, iflag)
      return medium, WF1, WF2, WFSize, FD
  #=================== MOR forward ============================#
    elseif morforward == true
      FD = Fdmtx(WFSize, medium, FDC, ext, iflag)
      return medium, WFSize, FD
    end
end

############ elastic layered medium with self-defined layer properties ############
function InitMod{T1,T2,T3,T4,T5,T6,T7<:Real}(vp::Array{T1,1},
  vs::Array{T2,1}, rho::Array{T3,1}, Depth::Array{T4,1},Horizon::T5, PeakFreq::T6, T::T7; ext = 10, iflag = 2, mortraining = false, morforward = false)

  medium = model(vp,vs,rho, Depth, Horizon, PeakFreq, T, ext, iflag)
  #%%%%%%%% wavefield initialization %%%%%%%%%
  WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
  WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
  WFSize = InitWF2DSize(medium.NDep, medium.NHor, ext, iflag)
  ############### Finite Difference Matrix ################
  FDC = FDCoeff(4)

  # =================== unsplit PML boundary ===================================#
    if mortraining == false && morforward == false
      #%%%%%%%% wavefield initialization %%%%%%%%%
      WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
      WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)

      FD = Fdmtx(WFSize, medium, FDC, ext)
      pml = DampBound(WFSize, medium, ext, iflag)
      return medium, WF1, WF2, WFSize, FD, pml
  # =================== split PML boundary ===================================#
    elseif mortraining == true && morforward == false
      WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
      WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)

      FD = Fdmtx(WFSize, medium, FDC, ext, iflag)
      return medium, WF1, WF2, WFSize, FD
  #=================== MOR forward ============================#
    elseif morforward == true
      FD = Fdmtx(WFSize, medium, FDC, ext, iflag)
      return medium, WFSize, FD
    end
end

############ medium properties input from Vp, Vs and density file ############
function InitMod{T1,T2,T3,T4<:Real}(vp_path::String,
  vs_path::String, rho_path::String,dz::T1, dx::T2, NDep::Int64, NHor::Int64,PeakFreq::T3, T::T4; ext = 10, iflag = 2, mortraining = false, morforward = false)

  medium = model(vp_path,vs_path,rho_path,dz,dx,NDep,NHor,PeakFreq,T,ext,iflag)
  #%%%%%%%% wavefield initialization %%%%%%%%%
  WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
  WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
  WFSize = InitWF2DSize(medium.NDep, medium.NHor, ext, iflag)
  ############### Finite Difference Matrix ################
  FDC = FDCoeff(4)

  # =================== unsplit PML boundary ===================================#
    if mortraining == false && morforward == false
      #%%%%%%%% wavefield initialization %%%%%%%%%
      WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
      WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)

      FD = Fdmtx(WFSize, medium, FDC, ext)
      pml = DampBound(WFSize, medium, ext, iflag)
      return medium, WF1, WF2, WFSize, FD, pml
  # =================== split PML boundary ===================================#
    elseif mortraining == true && morforward == false
      WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
      WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)

      FD = Fdmtx(WFSize, medium, FDC, ext, iflag)
      return medium, WF1, WF2, WFSize, FD
  #=================== MOR forward ============================#
    elseif morforward == true
      FD = Fdmtx(WFSize, medium, FDC, ext, iflag)
      return medium, WFSize, FD
    end
end

############ elastic medium with only Input Vp from file ############
function InitMod{T1,T2,T3,T4,T5<:Real}(vp_path::String,
  datatype::String, rho::T1, dz::T2, dx::T3, NDep::Int64, NHor::Int64, PeakFreq::T4, T::T5; ext = 10, iflag = 2, mortraining = false, morforward = false)

    medium = model(vp_path,datatype,rho,dz,dx,NDep,NHor,PeakFreq,T,ext,iflag)
    #%%%%%%%% wavefield initialization %%%%%%%%%
    WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
    WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
    WFSize = InitWF2DSize(medium.NDep, medium.NHor, ext, iflag)
    ############### Finite Difference Matrix ################
    FDC = FDCoeff(4)

    # =================== unsplit PML boundary ===================================#
      if mortraining == false && morforward == false
        #%%%%%%%% wavefield initialization %%%%%%%%%
        WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
        WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)

        FD = Fdmtx(WFSize, medium, FDC, ext)
        pml = DampBound(WFSize, medium, ext, iflag)
        return medium, WF1, WF2, WFSize, FD, pml
    # =================== split PML boundary ===================================#
      elseif mortraining == true && morforward == false
        WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
        WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)

        FD = Fdmtx(WFSize, medium, FDC, ext, iflag)
        return medium, WF1, WF2, WFSize, FD
    #=================== MOR forward ============================#
      elseif morforward == true
        FD = Fdmtx(WFSize, medium, FDC, ext, iflag)
        return medium, WFSize, FD
      end
end
#
# # =================== split PML boundary ===================================#
# ############ homogeneous acoustic medium ############
# function InitSLMod{T1,T2,T3,T4,T5,T6<:Real}(v_p::T1,
#   rho::T2,Depth::T3,Horizon::T4,PeakFreq::T5,T::T6; ext = 10, iflag = 2)
#
#   medium = model(v_p, rho, Depth, Horizon, PeakFreq, T, ext, iflag)
#   #%%%%%%%% wavefield initialization %%%%%%%%%
#   WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
#   WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
#   WFSize = InitWF2DSize(medium.NDep, medium.NHor, ext, iflag)
#   ############### Finite Difference Matrix ################
#   FDC = FDCoeff(4)
#
#   FD = SLFdmtx(WFSize, medium, FDC, ext, iflag)
#
#   return medium, WF1, WF2, WFSize, FD
# end
#
# ############ homogeneous elastic medium ############
# function InitSLMod{T1,T2,T3,T4,T5,T6,T7<:Real}(v_p::T1,
#   v_s::T2, rho::T3, Depth::T4, Horizon::T5, PeakFreq::T6, T::T7; ext = 10, iflag = 2)
#
#   medium = model(v_p, v_s, rho, Depth, Horizon, PeakFreq, T, ext, iflag)
#   #%%%%%%%% wavefield initialization %%%%%%%%%
#   WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
#   WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
#   WFSize = InitWF2DSize(medium.NDep, medium.NHor, ext, iflag)
#   ############### Finite Difference Matrix ################
#   FDC = FDCoeff(4)
#
#   FD = SLFdmtx(WFSize, medium, FDC, ext, iflag)
#
#   return medium, WF1, WF2, WFSize, FD
# end
#
# ############ elastic layered medium with constant density ############
# function InitSLMod{T1,T2,T3,T4,T5,T6<:Real}(Max_Min_v::Array{T1,1},
#   rho::T2, LayerNo::Int64,Depth::T3, Horizon::T4, PeakFreq::T5, T::T6; ext = 10, iflag = 2)
#
#   medium = model(Max_Min_v, rho, LayerNo, Depth, Horizon, PeakFreq, T, ext, iflag)
#   #%%%%%%%% wavefield initialization %%%%%%%%%
#   WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
#   WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
#   WFSize = InitWF2DSize(medium.NDep, medium.NHor, ext, iflag)
#   ############### Finite Difference Matrix ################
#   FDC = FDCoeff(4)
#
#   FD = SLFdmtx(WFSize, medium, FDC, ext, iflag)
#
#   return medium, WF1, WF2, WFSize, FD
# end
#
# ############ elastic layered medium with layerd density and constant layer depth ############
# function InitSLMod{T1,T2,T3,T4,T5,T6<:Real}(Max_Min_v::Array{T1,1},
#   Max_Min_rho::Array{T2,1}, LayerNo::Int64, Depth::T3, Horizon::T4, PeakFreq::T5, T::T6; ext = 10, iflag = 2)
#
#   medium = model(Max_Min_v, Max_Min_rho, LayerNo, Depth, Horizon, PeakFreq, T, ext, iflag)
#   #%%%%%%%% wavefield initialization %%%%%%%%%
#   WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
#   WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
#   WFSize = InitWF2DSize(medium.NDep, medium.NHor, ext, iflag)
#   ############### Finite Difference Matrix ################
#   FDC = FDCoeff(4)
#
#   FD = SLFdmtx(WFSize, medium, FDC, ext, iflag)
#
#   return medium, WF1, WF2, WFSize, FD
# end
#
# ############ elastic layered medium with self-defined layer properties ############
# function InitSLMod{T1,T2,T3,T4,T5,T6,T7<:Real}(vp::Array{T1,1},
#   vs::Array{T2,1}, rho::Array{T3,1}, Depth::Array{T4,1},Horizon::T5, PeakFreq::T6, T::T7; ext = 10, iflag = 2)
#
#   medium = model(vp,vs,rho, Depth, Horizon, PeakFreq, T, ext, iflag)
#   #%%%%%%%% wavefield initialization %%%%%%%%%
#   WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
#   WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
#   WFSize = InitWF2DSize(medium.NDep, medium.NHor, ext, iflag)
#   ############### Finite Difference Matrix ################
#   FDC = FDCoeff(4)
#
#   FD = SLFdmtx(WFSize, medium, FDC, ext, iflag)
#
#   return medium, WF1, WF2, WFSize, FD
# end
#
# ############ medium properties input from Vp, Vs and density file ############
# function InitSLMod{T1,T2,T3,T4<:Real}(vp_path::String,
#   vs_path::String, rho_path::String,dz::T1, dx::T2, NDep::Int64, NHor::Int64,PeakFreq::T3, T::T4; ext = 10, iflag = 2)
#
#   medium = model(vp_path,vs_path,rho_path,dz,dx,NDep,NHor,PeakFreq,T,ext,iflag)
#   #%%%%%%%% wavefield initialization %%%%%%%%%
#   WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
#   WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
#   WFSize = InitWF2DSize(medium.NDep, medium.NHor, ext, iflag)
#   ############### Finite Difference Matrix ################
#   FDC = FDCoeff(4)
#
#   FD = SLFdmtx(WFSize, medium, FDC, ext, iflag)
#
#   return medium, WF1, WF2, WFSize, FD
# end
#
# ############ elastic medium with only Input Vp from file ############
# function InitSLMod{T1,T2,T3,T4,T5<:Real}(vp_path::String,
#   datatype::String, rho::T1, dz::T2, dx::T3, NDep::Int64, NHor::Int64, PeakFreq::T4, T::T5; ext = 10, iflag = 2)
#
#     medium = model(vp_path,datatype,rho,dz,dx,NDep,NHor,PeakFreq,T,ext,iflag)
#     #%%%%%%%% wavefield initialization %%%%%%%%%
#     WF1 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
#     WF2 = InitWF2D(medium.NDep, medium.NHor, ext, iflag)
#     WFSize = InitWF2DSize(medium.NDep, medium.NHor, ext, iflag)
#     ############### Finite Difference Matrix ################
#     FDC = FDCoeff(4)
#
#     FD = SLFdmtx(WFSize, medium, FDC, ext, iflag)
#
#     return medium, WF1, WF2, WFSize, FD
# end
