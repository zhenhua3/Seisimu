#=== elastic ===#
#=== write MSB ===#
function wMSB(
    model::Union{elmod2d},
    sou::Array{source},
    sourcenumber::Int64,
    ParbinPath::Union{String, Void})

    nwf = model.nwf
    pkf = model.medium.pkf
    iflag = model.medium.iflag
    medium = model.medium
    ext = model.medium.ext

    if ParbinPath != nothing
        if ParbinPath[end-2:end] != "bin"
            error("Check output file extensions. It has to be in 'bin' format.")
        else
            fid = open(ParbinPath,"a+")
            write(fid,iflag)
            write(fid,ext)
            write(fid,nwf.BDntpp)
            write(fid,medium.pvel)
            write(fid,medium.svel)
            write(fid,medium.rho)
            write(fid,Float64(medium.dx))
            write(fid,Float64(medium.dz))
            write(fid,Float64(medium.dt))
            write(fid,medium.nT)
            write(fid,Float64(pkf))

            write(fid,sourcenumber)
            for i = 1:sourcenumber
                write(fid,sou[i].nloc)
            end
            for i = 1:sourcenumber
                write(fid,sou[i].waveform)
            end
            for i = 1:sourcenumber
                write(fid,sou[i].not)
            end
            for i = 1:sourcenumber
                if sou[i].tp == "expl"
                    write(fid,1)
                elseif sou[i].tp == "sfx"
                    write(fid,2)
                elseif sou[i].tp == "sfz"
                    write(fid,4)
                elseif sou[i].tp == "scx"
                    write(fid,5)
                elseif sou[i].tp == "scz"
                    write(fid,7)
                elseif sou[i].tp == "DC"
                    write(fid,8)
                elseif sou[i].tp == "wy"
                    write(fid,10)
                end
            end
            close(fid)
        end
    end
end



















#=== acoustic ===#
#=== write MSB ===#
function wMSB(
    model::Union{acmod2d},
    sou::Array{source},
    sourcenumber::Int64,
    ParbinPath::Union{String, Void})

    nwf = model.nwf
    pkf = model.medium.pkf
    iflag = model.medium.iflag
    medium = model.medium
    ext = model.medium.ext

    if ParbinPath != nothing
        if ParbinPath[end-2:end] != "bin"
            error("Check output file extensions. It has to be in 'bin' format.")
        else
            fid = open(ParbinPath,"a+")
            write(fid,iflag)
            write(fid,ext)
            write(fid,nwf.BDntpp)
            write(fid,medium.pvel)
            write(fid,medium.rho)
            write(fid,Float64(medium.dx))
            write(fid,Float64(medium.dz))
            write(fid,Float64(medium.dt))
            write(fid,medium.nT)
            write(fid,Float64(pkf))

            write(fid,sourcenumber)
            for i = 1:sourcenumber
                write(fid,sou[i].nloc)
            end
            for i = 1:sourcenumber
                write(fid,sou[i].waveform)
            end
            for i = 1:sourcenumber
                write(fid,sou[i].not)
            end
            for i = 1:sourcenumber
                if sou[i].tp == "expl"
                    write(fid,1)
                elseif sou[i].tp == "sfx"
                    write(fid,2)
                elseif sou[i].tp == "sfz"
                    write(fid,4)
                elseif sou[i].tp == "scx"
                    write(fid,5)
                elseif sou[i].tp == "scz"
                    write(fid,7)
                end
            end
            close(fid)
        end
    end
end
