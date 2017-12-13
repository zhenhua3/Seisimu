            #====================#
            #=== Medium State ===#
            #====================#

#=== elastic 2d ===#
 type elastic2d{T1,T2,T3,T4<:Real}
    pvel::Array{Float64,2}
    svel::Array{Float64,2}
    rho::Array{Float64,2}
    lambda::Array{Float64,2}
    mu::Array{Float64,2}
    dx::Float64
    dz::Float64
    dt::Float64
    DZ::T1
    HX::T2
    nDZ::Int64
    nHX::Int64
    BDnDZ::Int64
    BDnHX::Int64
    T::T3 # Total Simulation Time
    nT::Int64 # Total discretized sampling time
    pkf::T4
    ext::Int64
    iflag::Int64
end
#=== acoustic 2d ===#
 type acoustic2d{T1,T2,T3,T4<:Real}
    pvel::Array{Float64,2}
    rho::Array{Float64,2}
    lambda::Array{Float64,2}
    dx::Float64
    dz::Float64
    dt::Float64
    DZ::T1
    HX::T2
    nDZ::Int64
    nHX::Int64
    BDnDZ::Int64
    BDnHX::Int64
    T::T3 # Total Simulation Time
    nT::Int64 # Total discretized sampling time
    pkf::T4
    ext::Int64
    iflag::Int64
end

#=== elastic 3d ===#
 type elastic3d{T1,T2,T3,T4,T5<:Real}
    pvel::Array{Float32,3}
    svel::Array{Float32,3}
    rho::Array{Float32,3}
    lambda::Array{Float32,3}
    mu::Array{Float32,3}
    dx::Float32
    dy::Float32
    dz::Float32
    dt::Float32
    DZ::T1
    HX::T2
    HY::T3
    nDZ::Int64
    nHX::Int64
    nHY::Int64
    BDnDZ::Int64
    BDnHX::Int64
    BDnHY::Int64
    T::T4 # Total Simulation Time
    nT::Int64 # Total discretized sampling time
    pkf::T5
    ext::Int64
    iflag::Int64
end

#=== acoustic 3d ===#
 type acoustic3d{T1,T2,T3,T4,T5<:Real}
    pvel::Array{Float32,3}
    rho::Array{Float32,3}
    lambda::Array{Float32,3}
    dx::Float32
    dy::Float32
    dz::Float32
    dt::Float32
    DZ::T1
    HX::T2
    HY::T3
    nDZ::Int64
    nHX::Int64
    nHY::Int64
    BDnDZ::Int64
    BDnHX::Int64
    BDnHY::Int64
    T::T4 # Total Simulation Time
    nT::Int64 # Total discretized sampling time
    pkf::T5
    ext::Int64
    iflag::Int64
end

            #==========================#
            #=== Obsorbing boundary ===#
            #==========================#

#=== elastic 2d ===#
 type elbd2d
    bfull::Array{Float64,1}
    bhalf::Array{Float64,1}
    afull::Array{Float64,1}
    ahalf::Array{Float64,1}
    PVxBTxx::Array{Float64,2}
    PVxBTxz::Array{Float64,2}
    PVzBTzz::Array{Float64,2}
    PVzBTxz::Array{Float64,2}
    PTxxBVx::Array{Float64,2}
    PTxxBVz::Array{Float64,2}
    PTzzBVx::Array{Float64,2}
    PTzzBVz::Array{Float64,2}
    PTxzBVx::Array{Float64,2}
    PTxzBVz::Array{Float64,2}
end
#=== acoustic 2d ===#
 type acbd2d
    bfull::Array{Float64,1}
    bhalf::Array{Float64,1}
    afull::Array{Float64,1}
    ahalf::Array{Float64,1}
    PVxBTxx::Array{Float64,2}
    PVzBTzz::Array{Float64,2}
    PTxxBVx::Array{Float64,2}
    PTxxBVz::Array{Float64,2}
    PTzzBVx::Array{Float64,2}
    PTzzBVz::Array{Float64,2}
end
#=== elastic 3d ===#
 type elbd3d
    bfull::Array{Float32,1}
    bhalf::Array{Float32,1}
    afull::Array{Float32,1}
    ahalf::Array{Float32,1}
    PVzBTzz::Array{Float32,3}
    PVzBTxz::Array{Float32,3}
    PVzBTyz::Array{Float32,3}
    PVxBTxx::Array{Float32,3}
    PVxBTxz::Array{Float32,3}
    PVxBTxy::Array{Float32,3}
    PVyBTyy::Array{Float32,3}
    PVyBTyz::Array{Float32,3}
    PVyBTxy::Array{Float32,3}
    PTzzBVz::Array{Float32,3}
    PTzzBVx::Array{Float32,3}
    PTzzBVy::Array{Float32,3}
    PTxxBVz::Array{Float32,3}
    PTxxBVx::Array{Float32,3}
    PTxxBVy::Array{Float32,3}
    PTyyBVz::Array{Float32,3}
    PTyyBVx::Array{Float32,3}
    PTyyBVy::Array{Float32,3}
    PTxzBVx::Array{Float32,3}
    PTxzBVz::Array{Float32,3}
    PTyzBVz::Array{Float32,3}
    PTyzBVy::Array{Float32,3}
    PTxyBVx::Array{Float32,3}
    PTxyBVy::Array{Float32,3}
end
#=== ascoustic 3d ===#
 type acbd3d
    bfull::Array{Float32,1}
    bhalf::Array{Float32,1}
    afull::Array{Float32,1}
    ahalf::Array{Float32,1}
    PVzBTzz::Array{Float32,3}
    PVxBTxx::Array{Float32,3}
    PVyBTyy::Array{Float32,3}
    PTzzBVz::Array{Float32,3}
    PTzzBVx::Array{Float32,3}
    PTzzBVy::Array{Float32,3}
    PTxxBVz::Array{Float32,3}
    PTxxBVx::Array{Float32,3}
    PTxxBVy::Array{Float32,3}
    PTyyBVz::Array{Float32,3}
    PTyyBVx::Array{Float32,3}
    PTyyBVy::Array{Float32,3}
end


            #========================#
            #=== wavefield define ===#
            #========================#

#=== non split elastic wavefield 2d ===#
 type elwf2d
    txx :: Array{Float64,2}
    tzz :: Array{Float64,2}
    txz :: Array{Float64,2}
    vx :: Array{Float64,2}
    vz :: Array{Float64,2}
    tmp_tpp :: Array{Float64,2}
    tmp_txz :: Array{Float64,2}
    tmp_vx :: Array{Float64,2}
    tmp_vz :: Array{Float64,2}
end

#=== non split acoustic wavefield 2d ===#
 type acwf2d
    txx :: Array{Float64,2}
    tzz :: Array{Float64,2}
    vx :: Array{Float64,2}
    vz :: Array{Float64,2}
    tmp_tpp :: Array{Float64,2}
    tmp_vx :: Array{Float64,2}
    tmp_vz :: Array{Float64,2}
end

#=== non split elastic wavefield 3d ===#
 type elwf3d
  vz :: Array{Float32,3}
  vx :: Array{Float32,3}
  vy :: Array{Float32,3}
  tzz :: Array{Float32,3}
  txx :: Array{Float32,3}
  tyy :: Array{Float32,3}
  txz :: Array{Float32,3}
  txy :: Array{Float32,3}
  tyz :: Array{Float32,3}
  tmp_vz :: Array{Float32,2}
  tmp_vx :: Array{Float32,2}
  tmp_vy :: Array{Float32,2}
  tmp_tpp :: Array{Float32,2}
  tmp_txz :: Array{Float32,2}
  tmp_txy :: Array{Float32,2}
  tmp_tyz :: Array{Float32,2}
end

#=== non split acoustic wavefield 3d ===#
 type acwf3d
  vz :: Array{Float32,3}
  vx :: Array{Float32,3}
  vy :: Array{Float32,3}
  tzz :: Array{Float32,3}
  txx :: Array{Float32,3}
  tyy :: Array{Float32,3}
  tmp_vz :: Array{Float32,2}
  tmp_vx :: Array{Float32,2}
  tmp_vy :: Array{Float32,2}
  tmp_tpp :: Array{Float32,2}
end



            #======================#
            #=== Wavefield Size ===#
            #======================#

#=== elastic 2d ===#
 type nelwf2d
  BDntpp :: Array{Int64,1}
  BDntxz :: Array{Int64,1}
  BDnvx :: Array{Int64,1}
  BDnvz :: Array{Int64,1}
  ntpp :: Array{Int64,1}
  ntxz :: Array{Int64,1}
  nvx :: Array{Int64,1}
  nvz :: Array{Int64,1}
end
#=== acoustic 2d ===#
 type nacwf2d
  BDntpp :: Array{Int64,1}
  BDnvx :: Array{Int64,1}
  BDnvz :: Array{Int64,1}
  ntpp :: Array{Int64,1}
  nvx :: Array{Int64,1}
  nvz :: Array{Int64,1}
end
#=== elastic 3d ===#
 type nelwf3d
    BDnvz :: Array{Int64,1}
    BDnvx :: Array{Int64,1}
    BDnvy :: Array{Int64,1}
    BDntpp :: Array{Int64,1}
    BDntxz :: Array{Int64,1}
    BDntxy :: Array{Int64,1}
    BDntyz :: Array{Int64,1}
    nvz :: Array{Int64,1}
    nvx :: Array{Int64,1}
    nvy :: Array{Int64,1}
    ntpp :: Array{Int64,1}
    ntxz :: Array{Int64,1}
    ntxy :: Array{Int64,1}
    ntyz :: Array{Int64,1}
end
#=== acoustic 3d ===#
 type nacwf3d
    BDnvz :: Array{Int64,1}
    BDnvx :: Array{Int64,1}
    BDnvy :: Array{Int64,1}
    BDntpp :: Array{Int64,1}
    nvz :: Array{Int64,1}
    nvx :: Array{Int64,1}
    nvy :: Array{Int64,1}
    ntpp :: Array{Int64,1}
end

            #========================#
            #=== Wavefield Paras  ===#
            #========================#

#=== elastic 2d ===#
type elwfpara2d
   rho_vz::Array{Float64}
   rho_vx::Array{Float64}
   lamu::Array{Float64}
   mu_txz::Array{Float64}
end

#=== acoustic 2d ===#
type acwfpara2d
   rho_vz::Array{Float64}
   rho_vx::Array{Float64}
end

#=== elastic 3d ===#
type elwfpara3d
   rho_vz::Array{Float32}
   rho_vx::Array{Float32}
   rho_vy::Array{Float32}
   lamu::Array{Float32}
   mu_txz::Array{Float32}
   mu_txy::Array{Float32}
   mu_tyz::Array{Float32}
end

#=== acoustic 3d ===#
type acwfpara3d
   rho_vz::Array{Float32}
   rho_vx::Array{Float32}
   rho_vy::Array{Float32}
end

            #========================#
            #=== Simulation Model ===#
            #========================#


#=== non split elastic 2d ===#
type elmod2d
    medium::elastic2d
    wf::elwf2d
    nwf::nelwf2d
    fdc::Array{Float64,1}
    pml::elbd2d
    calpara::elwfpara2d
end

#=== non split acoustic 2d ===#
type acmod2d
    medium::acoustic2d
    wf::acwf2d
    nwf::nacwf2d
    fdc::Array{Float64,1}
    pml::acbd2d
    calpara::acwfpara2d
end

#=== non split elastic 3d ===#
type elmod3d
    medium::elastic3d
    wf::elwf3d
    nwf::nelwf3d
    fdc::Array{Float32,1}
    pml::elbd3d
    calpara::elwfpara3d
end

#=== non split acoustic 3d ===#
type acmod3d
    medium::acoustic3d
    wf::acwf3d
    nwf::nacwf3d
    fdc::Array{Float32,1}
    pml::acbd3d
    calpara::acwfpara3d
end
