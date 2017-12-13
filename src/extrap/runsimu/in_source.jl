#=== 2d ===#
function addsou!(wf::Union{elwf2d,acwf2d}, sou::Array{source}, it::Int64)
    sn = length(sou)
    # it : total time samplings
    # sn : sou number
        for isn = 1:sn
            wf.vx[sou[isn].BDnloc[1],sou[isn].BDnloc[2]] = wf.vx[sou[isn].BDnloc[1],sou[isn].BDnloc[2]] + (sou[isn].mt[1] + sou[isn].mt[3])*sou[isn].waveform[it]
            wf.vx[sou[isn].BDnloc[1],sou[isn].BDnloc[2]-1] = wf.vx[sou[isn].BDnloc[1],sou[isn].BDnloc[2]-1] - sou[isn].mt[1]*sou[isn].waveform[it]
            wf.vx[sou[isn].BDnloc[1]-1,sou[isn].BDnloc[2]] = wf.vx[sou[isn].BDnloc[1]-1,sou[isn].BDnloc[2]] - sou[isn].mt[3]*sou[isn].waveform[it]
            wf.vz[sou[isn].BDnloc[1],sou[isn].BDnloc[2]] = wf.vz[sou[isn].BDnloc[1],sou[isn].BDnloc[2]] + (sou[isn].mt[2] + sou[isn].mt[3])*sou[isn].waveform[it]
            wf.vz[sou[isn].BDnloc[1]-1,sou[isn].BDnloc[2]] = wf.vz[sou[isn].BDnloc[1]-1,sou[isn].BDnloc[2]] - sou[isn].mt[2]*sou[isn].waveform[it]
            wf.vz[sou[isn].BDnloc[1],sou[isn].BDnloc[2]-1] = wf.vz[sou[isn].BDnloc[1],sou[isn].BDnloc[2]-1] + sou[isn].mt[3]*sou[isn].waveform[it]
        end
end


#=== 3d ===#
function addsou!(wf::Union{elwf3d,acwf3d}, sou::Array{source}, it::Int64)
    sn = length(sou)
    # it : total time samplings
    # sn : sou number
        for isn = 1:sn
            wf.vx[sou[isn].BDnloc[1],sou[isn].BDnloc[2],sou[isn].BDnloc[3]] = wf.vx[sou[isn].BDnloc[1],sou[isn].BDnloc[2],sou[isn].BDnloc[3]] + (sou[isn].mt[1] + sou[isn].mt[4] + sou[isn].mt[5])*sou[isn].waveform[it]
            wf.vx[sou[isn].BDnloc[1],sou[isn].BDnloc[2]-1,sou[isn].BDnloc[3]] = wf.vx[sou[isn].BDnloc[1],sou[isn].BDnloc[2]-1,sou[isn].BDnloc[3]] - sou[isn].mt[1]*sou[isn].waveform[it]
            wf.vx[sou[isn].BDnloc[1],sou[isn].BDnloc[2],sou[isn].BDnloc[3]-1] = wf.vx[sou[isn].BDnloc[1],sou[isn].BDnloc[2],sou[isn].BDnloc[3]-1] - sou[isn].mt[4]*sou[isn].waveform[it]
            wf.vx[sou[isn].BDnloc[1]-1,sou[isn].BDnloc[2],sou[isn].BDnloc[3]] = wf.vx[sou[isn].BDnloc[1]-1,sou[isn].BDnloc[2],sou[isn].BDnloc[3]] - sou[isn].mt[5]*sou[isn].waveform[it]

            wf.vy[sou[isn].BDnloc[1],sou[isn].BDnloc[2],sou[isn].BDnloc[3]] = wf.vy[sou[isn].BDnloc[1],sou[isn].BDnloc[2],sou[isn].BDnloc[3]] + (sou[isn].mt[2] + sou[isn].mt[4] + sou[isn].mt[6])*sou[isn].waveform[it]
            wf.vy[sou[isn].BDnloc[1],sou[isn].BDnloc[2],sou[isn].BDnloc[3]-1] = wf.vy[sou[isn].BDnloc[1],sou[isn].BDnloc[2],sou[isn].BDnloc[3]-1] - sou[isn].mt[2]*sou[isn].waveform[it]
            wf.vy[sou[isn].BDnloc[1],sou[isn].BDnloc[2]-1,sou[isn].BDnloc[3]] = wf.vy[sou[isn].BDnloc[1],sou[isn].BDnloc[2]-1,sou[isn].BDnloc[3]] - sou[isn].mt[4]*sou[isn].waveform[it]
            wf.vy[sou[isn].BDnloc[1]-1,sou[isn].BDnloc[2],sou[isn].BDnloc[3]] = wf.vy[sou[isn].BDnloc[1]-1,sou[isn].BDnloc[2],sou[isn].BDnloc[3]] - sou[isn].mt[6]*sou[isn].waveform[it]


            wf.vz[sou[isn].BDnloc[1],sou[isn].BDnloc[2],sou[isn].BDnloc[3]] = wf.vz[sou[isn].BDnloc[1],sou[isn].BDnloc[2],sou[isn].BDnloc[3]] + (sou[isn].mt[3] + sou[isn].mt[5] + sou[isn].mt[6])*sou[isn].waveform[it]
            wf.vz[sou[isn].BDnloc[1]-1,sou[isn].BDnloc[2],sou[isn].BDnloc[3]] = wf.vz[sou[isn].BDnloc[1]-1,sou[isn].BDnloc[2],sou[isn].BDnloc[3]] - sou[isn].mt[2]*sou[isn].waveform[it]
            wf.vz[sou[isn].BDnloc[1],sou[isn].BDnloc[2]-1,sou[isn].BDnloc[3]] = wf.vz[sou[isn].BDnloc[1],sou[isn].BDnloc[2]-1,sou[isn].BDnloc[3]] - sou[isn].mt[5]*sou[isn].waveform[it]
            wf.vz[sou[isn].BDnloc[1],sou[isn].BDnloc[2],sou[isn].BDnloc[3]-1] = wf.vz[sou[isn].BDnloc[1],sou[isn].BDnloc[2],sou[isn].BDnloc[3]-1] - sou[isn].mt[6]*sou[isn].waveform[it]
        end
end
