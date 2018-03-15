# Check Output Component
function cpntlib2D()
    CpntLib = ("Vz","Vx","Tzz","Txx","P","Txz","Wy");
    println("Output Component Library: $CpntLib")
end

function CheckCpnt2D(
    rec::Union{Void,receiver},
    Slices::Union{Void,UnitRange{Int64},StepRange{Int64,Int64},Array{Int64,1},Int64},
    OptCpnt::Union{Void,Array{String}})
    if !(rec == nothing && Slices == nothing)
        if OptCpnt == nothing
            error("Please provide output component")
        else
            #CpntLib : Component Library
            CpntLib = ("Vz","Vx","Tzz","Txx","P","Txz","Wy","V","Tii");

            nOptCpnt = length(OptCpnt);
            # Check if Output component is in Component Library
            bool_rslt = map(x->x in CpntLib, OptCpnt);
            num_rslt = zeros(nOptCpnt);
            for i in 1:nOptCpnt
                num_rslt[i] = bool_rslt[i]*1
            end
            if countnz(num_rslt)!=nOptCpnt
                for i in 1:nOptCpnt
                    if bool_rslt[i] == false
                        println("$(OptCpnt[i]) is not a valid output component.")
                    end
                end
                error("Check output component library by typing 'cpntlib2D()'.")
            else println("Output Components Check Pass!")
            end
        end
    end
end

function message(
    rec::Union{Void,receiver},
    RecDataPath::Union{String,Void},
    Slices::Union{Void,UnitRange{Int64},StepRange{Int64,Int64},Array{Int64,1},Int64},
    WFDataPath::Union{String,Void},
    OptCpnt::Union{Void,Array{String}})

    if rec == nothing
        println("Warning: No recordings are output")
    else
        if RecDataPath == nothing
            error("Please provide recording output path")
        else
            println("Receiver locations are successfully initiated")
        end
    end
    if Slices == nothing
        println("Warning: No wavefields are output")
    else
        if WFDataPath == nothing
            error("Please provide wavefields output path")
        else
            println("Wavefield slices are successfully initiated")
        end
    end
    CheckCpnt2D(rec,Slices,OptCpnt)
end


function message(
    Slices::Union{Void,UnitRange{Int64},StepRange{Int64,Int64},Array{Int64,1},Int64},
    WFDataPath::Union{String,Void})

    if Slices == nothing
        error("Please give time slices for model order reduction")
    else
        if WFDataPath == nothing
            error("Please provide wavefields output path")
        else
            println("Wavefield slices are successfully initiated")
        end
    end
end
