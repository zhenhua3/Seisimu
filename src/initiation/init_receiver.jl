type receiver
    nr::Int64
    loc::Array{<:Real}
    nloc::Array{Int64}
    BDnloc::Array{Int64}
end

#=== receiver ===#
function initrec(
    nr::Int64,
    loc::Array{<:Real},
    model::Union{elmod2d,acmod2d})

    ext = model.medium.ext
    iflag = model.medium.iflag
    if length(loc) == 2 #=== receiver 2d ===#
        dz = model.medium.dz
        dx = model.medium.dx
        nloc = zeros(Int64,nr,2)
        for i = 1:nr
            nloc[i,1] = Int64(round(loc[i,1]./dz))
            nloc[i,2] = Int64(round(loc[i,2]./dx))
        end
        if iflag == 1 # free surface
            BDnloc = [nloc[:,1] nloc[:,2].+ext]
        elseif iflag == 2 # unlimited space
            BDnloc = [nloc[:,1].+ext nloc[:,2].+ext]
        end
    elseif length(loc) == 3 #=== receiver 3d ===#
        dz = model.medium.dz
        dx = model.medium.dx
        dy = model.medium.dy
        nloc = zeros(Int64,nr,3)
        for i = 1:nr
            nloc[i,1] = Int64(round(loc[i,1]./dz))
            nloc[i,2] = Int64(round(loc[i,2]./dx))
            nloc[i,3] = Int64(round(loc[i,3]./dy))
        end
        if iflag == 1 # free surface
            BDnloc = [nloc[:,1] nloc[:,2].+ext nloc[:,3].+ext]
        elseif iflag == 2 # unlimited space
            BDnloc = [nloc[:,1].+ext nloc[:,2].+ext nloc[:,3].+ext]
        end
    end
    return receiver(nr,loc,nloc,BDnloc)
end
