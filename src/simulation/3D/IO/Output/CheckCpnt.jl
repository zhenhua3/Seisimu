# Check Output Component
function cpntlib()
    CpntLib = ("Vz","Vx","Vy","V",
          "Tzz","Txx","Tyy","Tii","P",
          "Txz","Tyz","Txy","Tij",
          "Wz","Wx","Wy","W","Full");
    println("Output Component Library: $CpntLib")
end

function CheckCpnt(OptCpnt::Array{String})
    #CpntLib : Component Library
    CpntLib = ("Vz","Vx","Vy","V",
          "Tzz","Txx","Tyy","Tii","P",
          "Txz","Tyz","Txy","Tij",
          "Wz","Wx","Wy","W","Full");

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
        error("Check output component library by typing 'cpntlib()'.")
    else println("Output Components Check Pass!")
    end
end
