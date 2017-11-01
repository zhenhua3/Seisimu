type receiver
    rn::Int64
    loc::Array{<:Real}
    nloc::Array{Int64}
    BDnloc::Array{Int64}
end

function initrec(
    rn::Int64,
    loc::Array{<:Real},
    OptCpnt::Array{String},
    model::Union{spmod2d,nspmod2d})

    dz = model.medium.dz
    dx = model.medium.dx
    ext = model.medium.ext
    iflag = model.medium.iflag

    nloc = zeros(Int64,rn,2)
    for i = 1:rn
        nloc[i,1] = Int64(round(loc[i,1]./dz))
        nloc[i,2] = Int64(round(loc[i,2]./dx))
    end
    if iflag == 1 # free surface
        BDnloc = [nloc[:,1] nloc[:,2].+ext]
    elseif iflag == 2 # unlimited space
        BDnloc = [nloc[:,1].+ext nloc[:,2].+ext]
    end
    return receiver(rn,loc,nloc,BDnloc)
end
