type Model{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11<:Real}
    VP::Array{T1,2}
    VS::Array{T2,2}
    Rho::Array{T3,2}
    Lambda::Array{T4,2}
    Mu::Array{T5,2}
    dx::T6
    dz::T7
    dt::T8
    Dep::T9
    Hor::T10
    NDep::Int64
    NHor::Int64
    PML_NDep::Int64
    PML_NHor::Int64
    TST::T11 # Total Simulation Time
    Tn::Int64 # Total discretized sampling time
end

# This function calculates medium parameters for self-defined model or Input model
# Output:
    # Meshed P model with PML boudary
    # Meshed S model with PML boudary
    # Meshed Rho model with PML boundary
    # Meshed Lambda model with PML boudary
    # Meshed Mu model with PML boundary
    # dx : x spatial interval
    # dz : z spatial interval
    # dt : time interval
    # NDep : spatial sampling number in depth
    # NHor : spatial sampling number in horizontal distance
    # PML_NDep : spatial sampling number in depth with PML extend
    # PML_NHor : spatial sampling number in horizontal distance with PML extend

function ModExpand(Par::AbstractArray, ext::Int64, iflag::Int64)
            (NDep,NHor) = size(Par)
    if iflag == 1 #free surface
                temp_vertical2 = repmat(Par[end:end,:],ext,1)
                Para = vcat(Par, temp_vertical2)
                temp_horizon1 = repmat(Para[:,1],1,ext)
                temp_horizon2 = repmat(Para[:,end],1,ext)
                Para = hcat(temp_horizon2,Para,temp_horizon1)
                PML_NDep = NDep + ext
                PML_NHor = NHor + 2*ext
        elseif iflag == 2  # unlimited medium
                temp_vertical1 = repmat(Par[1:1,:],ext,1)
                temp_vertical2 = repmat(Par[end:end,:],ext,1)
                Para = vcat(temp_vertical1,Par,temp_vertical2)
                temp_horizon1 = repmat(Para[:,1],1,ext)
                temp_horizon2 = repmat(Para[:,end],1,ext)
                Para = hcat(temp_horizon1,Para,temp_horizon2)
                PML_NDep = NDep + 2*ext
                PML_NHor = NHor + 2*ext
            end
        return Para, PML_NDep, PML_NHor
end

############ homogeneous acoustic medium ############
function model{T1,T2,T3,T4,T5,T6<:Real}(v_p::T1,rho::T2,Depth::T3,Horizon::T4,PeakFreq::T5,T::T6, ext::Int64, iflag::Int64)
# Input:
    # v_p : homogeneous P-wave velocity
    # rho : Rho
    # Depth : Total depth of the model
    # Horizon : Total horizontal distance of the model
    # PeakFreq : Peak frequency for the wavelet
    # medium : Output
    # T : Total simulation time
    # ext : PML boundary for oneside with default value 15
    # iflag : 1 : free surface; 2: unlimited surface default value 2

    # Spatial interval and time interval
    dz = v_p / 8 / PeakFreq
    dx = v_p / 8 / PeakFreq
    dt = 0.406 * dz / v_p
    TST= T
    Tn = Int64(round(TST/dt))

    # discretize model
    NDep = Int64(round(Depth/dz))
    NHor = Int64(round(Horizon/dx))
    Vp = zeros(NDep,NHor).+v_p
    Vs = zeros(NDep,NHor)
    Rho = zeros(NDep,NHor).+rho

    VP, PML_NDep, PML_NHor = ModExpand(Vp, ext, iflag) # 1 : free surface; 2: unlimited surface
    VS, PML_NDep, PML_NHor = ModExpand(Vs, ext, iflag)
    Rho, PML_NDep, PML_NHor = ModExpand(Rho, ext, iflag)
    Lambda = VP.^2.*Rho - 2*VS.^2.*Rho
    Mu = VS.^2.*Rho
    Dep = Depth
    Hor = Horizon
    medium = Model(VP,VS,Rho,Lambda,Mu,dx,dz,dt,Dep,Hor,NDep,NHor,PML_NDep,PML_NHor,TST,Tn)
    return medium
end

############ homogeneous elastic medium ############
function model{T1,T2,T3,T4,T5,T6,T7<:Real}(v_p::T1, v_s::T2, rho::T3, Depth::T4, Horizon::T5, PeakFreq::T6, T::T7, ext::Int64, iflag::Int64)
# Input:
    # v_p : homogeneous P-wave velocity
    # v_s : homogeneous S-wave velocity
    # rho : Rho
    # Depth : Total depth of the model
    # Horizon : Total horizontal distance of the model
    # PeakFreq : Peak frequency for the wavelet
    # medium : Output
    # T : Total simulation time
    # ext : PML boundary for oneside
    # iflag : 1 : free surface; 2: unlimited surface

    # Spatial interval and time interval
    if v_s == 0
      error("please remove vs from argument")
    end
    dz = v_s / 8 / PeakFreq
    dx = v_s / 8 / PeakFreq
    dt = 0.406 * dz / v_p
    TST= T
    Tn = Int64(round(TST/dt))

    # discretize model
    NDep = Int64(round(Depth/dz))
    NHor = Int64(round(Horizon/dx))
    Vp = zeros(NDep,NHor).+v_p
    Vs = zeros(NDep,NHor).+v_s
    Rho = zeros(NDep,NHor).+rho

    VP, PML_NDep, PML_NHor = ModExpand(Vp, ext, iflag) # 1 : free surface; 2: unlimited surface
    VS, PML_NDep, PML_NHor = ModExpand(Vs, ext, iflag)
    Rho, PML_NDep, PML_NHor = ModExpand(Rho, ext, iflag)
    Lambda = VP.^2.*Rho - 2*VS.^2.*Rho
    Mu = VS.^2.*Rho
    Dep = Depth
    Hor = Horizon
    medium = Model(VP,VS,Rho,Lambda,Mu,dx,dz,dt,Dep,Hor,NDep,NHor,PML_NDep,PML_NHor,TST,Tn)
    return medium
end

############ elastic layered medium with constant density ############
function model{T1,T2,T3,T4,T5,T6<:Real}(Max_Min_v::Array{T1,1}, rho::T2, LayerNo::Int64,Depth::T3, Horizon::T4, PeakFreq::T5, T::T6, ext::Int64, iflag::Int64)
# Input:
    # Max_Min_v : Maximum and Minimum P and S-wave velocity [Min_p, Max_p, Min_s, Max_s]
    # rho : constant Rho
    # LayerNo : Number of layers
    # Depth : Total depth of the model
    # Horizon : Total horizontal distance of the model
    # PeakFreq : Peak frequency for the wavelet
    # T : Total simulation time
    # medium : Output
    # ext : PML boundary for oneside
    # iflag : 1 : free surface; 2: unlimited surface

    if Max_Min_v[2]<Max_Min_v[1]
        error("Please check the order of your input P-wave velocity ")
    elseif Max_Min_v[4]<Max_Min_v[3]
        error("Please check the order of your input S-wave velocity ")
    elseif minimum(Max_Min_v) == 0
       error("velocity cannot be zero.")
    end

    max_vp = Max_Min_v[2]
    min_vp = Max_Min_v[1]
    max_vs = Max_Min_v[4]
    min_vs = Max_Min_v[3]
    min_v = minimum(Max_Min_v)

    # Spatial interval and time interval
    dz = min_v / 8 / PeakFreq
    dx = min_v / 8 / PeakFreq
    dt = 0.406 * dz / max_vp
    TST= T
    Tn = Int64(round(TST/dt))

    # discretize model
    NDep = Int64(round(Depth/dz))
    NHor = Int64(round(Horizon/dx))
    Vp = zeros(NDep,NHor)
    Vs = zeros(NDep,NHor)
    Rho = zeros(NDep,NHor).+rho

    # velocity interval and depth interval
    dvp = (max_vp - min_vp)/(LayerNo-1)
    dvs = (max_vs - min_vs)/(LayerNo-1)
    dDepth = Depth/LayerNo # depth of each layer

    vp = zeros(LayerNo)
    vs = zeros(LayerNo)

    for iL = 1 : LayerNo-1
    vp[iL] = (iL-1)*dvp + min_vp
    vs[iL] = (iL-1)*dvs + min_vs

    Vp[((iL-1)*Int64(round(dDepth/dz))+1):((iL-1)*Int64(round(dDepth/dz))+Int64(round(dDepth/dz))),:] = vp[iL]
    Vs[((iL-1)*Int64(round(dDepth/dz))+1):((iL-1)*Int64(round(dDepth/dz))+Int64(round(dDepth/dz))),:] = vs[iL]
    end
    iL = LayerNo
    vp[iL] = (iL-1)*dvp + min_vp
    vs[iL] = (iL-1)*dvs + min_vs

    Vp[((iL-1)*Int64(round(dDepth/dz))+1):NDep,:] = vp[iL]
    Vs[((iL-1)*Int64(round(dDepth/dz))+1):NDep,:] = vs[iL]

    VP, PML_NDep, PML_NHor = ModExpand(Vp, ext, iflag) # 1 : free surface; 2: unlimited surface
    VS, PML_NDep, PML_NHor = ModExpand(Vs, ext, iflag)
    Rho, PML_NDep, PML_NHor = ModExpand(Rho, ext, iflag)
    Lambda = VP.^2.*Rho - 2*VS.^2.*Rho
    Mu = VS.^2.*Rho
    Dep = Depth
    Hor = Horizon
    medium = Model(VP,VS,Rho,Lambda,Mu,dx,dz,dt,Dep,Hor,NDep,NHor,PML_NDep,PML_NHor,TST,Tn)
    return medium
end

############ elastic layered medium with layerd density and identical layer depth ############
function model{T1,T2,T3,T4,T5,T6<:Real}(Max_Min_v::Array{T1,1}, Max_Min_rho::Array{T2,1}, LayerNo::Int64, Depth::T3, Horizon::T4, PeakFreq::T5, T::T6, ext::Int64, iflag::Int64)
# Input:
    # Max_Min_v : Maximum and Minimum P and S-wave velocity [Min_p, Max_p, Min_s, Max_s]
    # Max_Min_rho : Maximum and Minimum Rho [Min_rho, Max_rho]
    # LayerNo : Number of layers
    # Depth : Total depth of the model
    # Horizon : Total horizontal distance of the model
    # PeakFreq : Peak frequency for the wavelet
    # T : Total simulation time
    # medium : Output
    # ext : PML boundary for oneside
    # iflag : 1 : free surface; 2: unlimited surface

    if Max_Min_v[2]<Max_Min_v[1]
        error("Please check the order of your input P-wave velocity ")
    elseif Max_Min_v[4]<Max_Min_v[3]
        error("Please check the order of your input S-wave velocity ")
    end

    max_vp = Max_Min_v[2]
    min_vp = Max_Min_v[1]
    max_vs = Max_Min_v[4]
    min_vs = Max_Min_v[3]
    min_v = minimum(Max_Min_v)
    max_rho = Max_Min_rho[1]
    min_rho = Max_Min_rho[2]

    # Spatial interval and time interval
    dz = min_v / 8 / PeakFreq
    dx = min_v / 8 / PeakFreq
    dt = 0.406 * dz / max_vp
    TST= T
    Tn = Int64(round(TST/dt))

    # discretize model
    NDep = Int64(round(Depth/dz))
    NHor = Int64(round(Horizon/dx))
    Vp = zeros(NDep,NHor)
    Vs = zeros(NDep,NHor)
    Rho = zeros(NDep,NHor)

    # velocity interval and depth interval
    dvp = (max_vp - min_vp)/(LayerNo-1)
    dvs = (max_vs - min_vs)/(LayerNo-1)
    drho = (max_rho - min_rho)/(LayerNo-1)
    dDepth = Depth/LayerNo # depth of each layer

    vp = zeros(LayerNo)
    vs = zeros(LayerNo)
    rho = zeros(LayerNo)

    for iL = 1 : LayerNo-1
    vp[iL] = (iL-1)*dvp + min_vp
    vs[iL] = (iL-1)*dvs + min_vs
    rho[iL] = (iL-1)*drho + min_rho

    Vp[((iL-1)*Int64(round(dDepth/dz))+1):((iL-1)*Int64(round(dDepth/dz))+Int64(round(dDepth/dz))),:] = vp[iL]
    Vs[((iL-1)*Int64(round(dDepth/dz))+1):((iL-1)*Int64(round(dDepth/dz))+Int64(round(dDepth/dz))),:] = vs[iL]
    Rho[((iL-1)*Int64(round(dDepth/dz))+1):((iL-1)*Int64(round(dDepth/dz))+Int64(round(dDepth/dz))),:] = rho[iL]
    end
    iL = LayerNo
    vp[iL] = (iL-1)*dvp + min_vp
    vs[iL] = (iL-1)*dvs + min_vs
    rho[iL] = (iL-1)*drho + min_rho

    Vp[((iL-1)*Int64(round(dDepth/dz))+1):NDep,:] = vp[iL]
    Vs[((iL-1)*Int64(round(dDepth/dz))+1):NDep,:] = vs[iL]
    Rho[((iL-1)*Int64(round(dDepth/dz))+1):NDep,:] = rho[iL]

    VP, PML_NDep, PML_NHor = ModExpand(Vp, ext, iflag) # 1 : free surface; 2: unlimited surface
    VS, PML_NDep, PML_NHor = ModExpand(Vs, ext, iflag)
    Rho, PML_NDep, PML_NHor = ModExpand(Rho, ext, iflag)
    Lambda = VP.^2.*Rho - 2*VS.^2.*Rho
    Mu = VS.^2.*Rho
    Dep = Depth
    Hor = Horizon
    medium = Model(VP,VS,Rho,Lambda,Mu,dx,dz,dt,Dep,Hor,NDep,NHor,PML_NDep,PML_NHor,TST,Tn)
    return medium
end

############ elastic layered medium with self-defined layer properties ############
function model{T1,T2,T3,T4,T5,T6,T7<:Real}(vp::Array{T1,1}, vs::Array{T2,1}, rho::Array{T3,1}, Depth::Array{T4,1},Horizon::T5, PeakFreq::T6, T::T7, ext::Int64, iflag::Int64)
# Input:
    # vp : p-wave velocity for each layer [vp1, vp2, vp3, ...]
    # vs : s-wave velocity for each layer [vs1, vs2, vs3, ...]
    # rho : rho for each layer [rho1, rho2, rho3, ...]
    # Depth : depth for each layer [depth1, depth2, depth3, ...]
    # Horizon : horzontal distance
    # PeakFreq : Peak frequency for the wavelet
    # T : Total simulation time
    # medium : Output
    # ext : PML boundary for oneside
    # iflag : 1 : free surface; 2: unlimited surface

    max_vp = maximum(vp)
    min_vp = minimum(vp)
    max_vs = maximum(vs)
    min_vs = minimum(vs)
    max_rho = maximum(rho)
    min_rho = minimum(rho)

    # Spatial interval and time interval
    dz = min_vs / 20 / PeakFreq
    dx = min_vs / 20 / PeakFreq
    dt = 0.406 * dz / max_vp
    TST= T
    Tn = Int64(round(TST/dt))

    # discretize model
    NDep = 0
  for i = 1:length(vp)
    NDep = NDep + Int64(round(Depth[i]/dz))
  end

    NHor = Int64(round(Horizon/dx))
    Vp = zeros(NDep,NHor)
    Vs = zeros(NDep,NHor)
    Rho = zeros(NDep,NHor)

    # velocity interval and depth interval
    LayerNo = length(vp)
    iDepth = 0

    for iL = 1 : LayerNo
    Vp[(iDepth+1):(iDepth+Int64(round(Depth[iL]/dz))),:] = vp[iL]
    Vs[(iDepth+1):(iDepth+Int64(round(Depth[iL]/dz))),:] = vs[iL]
    Rho[(iDepth+1):(iDepth+Int64(round(Depth[iL]/dz))),:] = rho[iL]
    iDepth = iDepth + Int64(round(Depth[iL]/dz))
    end

    VP, PML_NDep, PML_NHor = ModExpand(Vp, ext, iflag) # 1 : free surface; 2: unlimited surface
    VS, PML_NDep, PML_NHor = ModExpand(Vs, ext, iflag)
    Rho, PML_NDep, PML_NHor = ModExpand(Rho, ext, iflag)
    Lambda = VP.^2.*Rho - 2*VS.^2.*Rho
    Mu = VS.^2.*Rho
    Dep = sum(Depth)
    Hor = Horizon
    medium = Model(VP,VS,Rho,Lambda,Mu,dx,dz,dt,Dep,Hor,NDep,NHor,PML_NDep,PML_NHor,TST,Tn)
    return medium
end

############ medium properties input from Vp, Vs and density file ############
function model{T1,T2,T3,T4<:Real}(vp_path::String, vs_path::String, rho_path::String,dz::T1, dx::T2, NDep::Int64, NHor::Int64,PeakFreq::T3, T::T4, ext::Int64, iflag::Int64)
# Input:
    # vp_path : p-wave velocity file
    # vs_path : s-wave velocity file
    # rho_path : rho file
    # dz : spatial interval
    # dx : spatial interval
    # Nz : discretized number of depth
    # Nx : discretized number of horizon
    # PeakFreq : Peak frequency for the wavelet
    # T : Total simulation time
    # medium : Output
    # ext : PML boundary for oneside
    # iflag : 1 : free surface; 2: unlimited surface

    fid = open(vp_path,"r")
    vp = read(fid,Float64, NDep*NHor,1)
    close(fid)

    fid = open(vs_path,"r")
    vs = read(fid,Float64, NDep*NHor,1)
    close(fid)

    fid = open(rho_path,"r")
    rho = read(fid,Float64, NDep*NHor,1)
    close(fid)

    max_vp = maximum(vp)
    min_vp = minimum(vp)
    max_vs = maximum(vs)
    min_vs = minimum(vs)
    max_rho = maximum(rho)
    min_rho = minimum(rho)

    Vp = reshape(vp, NDep,NHor)
    Vs = reshape(vs, NDep,NHor)
    Rho = reshape(rho, NDep,NHor)

    # Spatial interval and time interval
    dz = dz
    dx = dx
    dt = 0.406 * dz / max_vp
    TST= T
    Tn = Int64(round(TST/dt))

    NDep = NDep
    NHor = NHor

    VP, PML_NDep, PML_NHor = ModExpand(Vp, ext, iflag) # 1 : free surface; 2: unlimited surface
    VS, PML_NDep, PML_NHor = ModExpand(Vs, ext, iflag)
    Rho, PML_NDep, PML_NHor = ModExpand(Rho, ext, iflag)
    Lambda = VP.^2.*Rho - 2*VS.^2.*Rho
    Mu = VS.^2.*Rho
    Dep = Depth
    Hor = Horizon
    medium = Model(VP,VS,Rho,Lambda,Mu,dx,dz,dt,Dep,Hor,NDep,NHor,PML_NDep,PML_NHor,TST,Tn)
    return medium
end

############  elastic medium with only Input Vp from file ############
function model{T1,T2,T3,T4,T5<:Real}(vp_path::String, datatype::String, rho::T1, dz::T2, dx::T3,NDep::Int64, NHor::Int64, PeakFreq::T4, T::T5, ext::Int64, iflag::Int64)
# Input:
    # vp_path : p-wave velocity file currently can read "bin" or "dat" file
    # rho_path : rho file
    # dz : spatial interval
    # dx : spatial interval
    # Nz : discretized number of depth
    # Nx : discretized number of horizon
    # PeakFreq : Peak frequency for the wavelet
    # T : Total simulation time
    # medium : Output
    # ext : PML boundary for oneside
    # iflag : 1 : free surface; 2: unlimited surface
if datatype == "bin"
    fid = open(vp_path,"r")
    vp = read(fid,Float64, NDep*NHor,1)
    close(fid)
elseif datatype == "dat"
  fid = open(vp_path,"r")
  vp = readdlm(fid)
  close(fid)
end

    vs = vp/sqrt(3)

    rho = zeros(NDep,NHor)+rho

    max_vp = maximum(vp)
    min_vp = minimum(vp)
    max_vs = maximum(vs)
    min_vs = minimum(vs)
    max_rho = maximum(rho)
    min_rho = minimum(rho)

    Vp = reshape(vp, NDep,NHor)
    Vs = reshape(vs, NDep,NHor)
    Rho = reshape(rho, NDep,NHor)

    # Spatial interval and time interval
    dz = dz
    dx = dx
    dt = 0.406 * dz / max_vp
    TST= T
    Tn = Int64(round(TST/dt))

    NDep = NDep
    NHor = NHor

    VP, PML_NDep, PML_NHor = ModExpand(Vp, ext, iflag) # 1 : free surface; 2: unlimited surface
    VS, PML_NDep, PML_NHor = ModExpand(Vs, ext, iflag)
    Rho, PML_NDep, PML_NHor = ModExpand(Rho, ext, iflag)
    Lambda = VP.^2.*Rho - 2*VS.^2.*Rho
    Mu = VS.^2.*Rho
    Dep = NDep*dz
    Hor = NHor*dx
    medium = Model(VP,VS,Rho,Lambda,Mu,dx,dz,dt,Dep,Hor,NDep,NHor,PML_NDep,PML_NHor,TST,Tn)
    return medium
end
