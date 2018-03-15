#=== 2d ===# # add moment tensor source #
function addmt!(wf::Union{elwf2d,acwf2d}, sou::MTsource, it::Int64)

    wf.vx[sou.BDnloc[1],sou.BDnloc[2]] = wf.vx[sou.BDnloc[1],sou.BDnloc[2]] + (sou.mt[1] + sou.mt[3])*sou.waveform[it]
    wf.vx[sou.BDnloc[1],sou.BDnloc[2]-1] = wf.vx[sou.BDnloc[1],sou.BDnloc[2]-1] - sou.mt[1]*sou.waveform[it]
    wf.vx[sou.BDnloc[1]-1,sou.BDnloc[2]] = wf.vx[sou.BDnloc[1]-1,sou.BDnloc[2]] - sou.mt[3]*sou.waveform[it]
    wf.vz[sou.BDnloc[1],sou.BDnloc[2]] = wf.vz[sou.BDnloc[1],sou.BDnloc[2]] + (sou.mt[2] + sou.mt[3])*sou.waveform[it]
    wf.vz[sou.BDnloc[1]-1,sou.BDnloc[2]] = wf.vz[sou.BDnloc[1]-1,sou.BDnloc[2]] - sou.mt[2]*sou.waveform[it]
    wf.vz[sou.BDnloc[1],sou.BDnloc[2]-1] = wf.vz[sou.BDnloc[1],sou.BDnloc[2]-1] + sou.mt[3]*sou.waveform[it]

end

#=== 2d ===# # add single force source #
function addsf!(wf::Union{elwf2d,acwf2d}, sou::SFsource, it::Int64)

    wf.vx[sou.BDnloc[1],sou.BDnloc[2]] = wf.vx[sou.BDnloc[1],sou.BDnloc[2]] + sou.coeff[1]*sou.waveform[it]
    wf.vz[sou.BDnloc[1],sou.BDnloc[2]] = wf.vz[sou.BDnloc[1],sou.BDnloc[2]] + sou.coeff[2]*sou.waveform[it]

end

#=== 3d ===#
function addmt!(wf::Union{elwf3d,acwf3d}, sou::MTsource, it::Int64)

    wf.vx[sou.BDnloc[1],sou.BDnloc[2],sou.BDnloc[3]] = wf.vx[sou.BDnloc[1],sou.BDnloc[2],sou.BDnloc[3]] + (sou.mt[1] + sou.mt[4] + sou.mt[5])*sou.waveform[it]
    wf.vx[sou.BDnloc[1],sou.BDnloc[2]-1,sou.BDnloc[3]] = wf.vx[sou.BDnloc[1],sou.BDnloc[2]-1,sou.BDnloc[3]] - sou.mt[1]*sou.waveform[it]
    wf.vx[sou.BDnloc[1],sou.BDnloc[2],sou.BDnloc[3]-1] = wf.vx[sou.BDnloc[1],sou.BDnloc[2],sou.BDnloc[3]-1] - sou.mt[4]*sou.waveform[it]
    wf.vx[sou.BDnloc[1]-1,sou.BDnloc[2],sou.BDnloc[3]] = wf.vx[sou.BDnloc[1]-1,sou.BDnloc[2],sou.BDnloc[3]] - sou.mt[5]*sou.waveform[it]

    wf.vy[sou.BDnloc[1],sou.BDnloc[2],sou.BDnloc[3]] = wf.vy[sou.BDnloc[1],sou.BDnloc[2],sou.BDnloc[3]] + (sou.mt[2] + sou.mt[4] + sou.mt[6])*sou.waveform[it]
    wf.vy[sou.BDnloc[1],sou.BDnloc[2],sou.BDnloc[3]-1] = wf.vy[sou.BDnloc[1],sou.BDnloc[2],sou.BDnloc[3]-1] - sou.mt[2]*sou.waveform[it]
    wf.vy[sou.BDnloc[1],sou.BDnloc[2]-1,sou.BDnloc[3]] = wf.vy[sou.BDnloc[1],sou.BDnloc[2]-1,sou.BDnloc[3]] - sou.mt[4]*sou.waveform[it]
    wf.vy[sou.BDnloc[1]-1,sou.BDnloc[2],sou.BDnloc[3]] = wf.vy[sou.BDnloc[1]-1,sou.BDnloc[2],sou.BDnloc[3]] - sou.mt[6]*sou.waveform[it]

    wf.vz[sou.BDnloc[1],sou.BDnloc[2],sou.BDnloc[3]] = wf.vz[sou.BDnloc[1],sou.BDnloc[2],sou.BDnloc[3]] + (sou.mt[3] + sou.mt[5] + sou.mt[6])*sou.waveform[it]
    wf.vz[sou.BDnloc[1]-1,sou.BDnloc[2],sou.BDnloc[3]] = wf.vz[sou.BDnloc[1]-1,sou.BDnloc[2],sou.BDnloc[3]] - sou.mt[2]*sou.waveform[it]
    wf.vz[sou.BDnloc[1],sou.BDnloc[2]-1,sou.BDnloc[3]] = wf.vz[sou.BDnloc[1],sou.BDnloc[2]-1,sou.BDnloc[3]] - sou.mt[5]*sou.waveform[it]
    wf.vz[sou.BDnloc[1],sou.BDnloc[2],sou.BDnloc[3]-1] = wf.vz[sou.BDnloc[1],sou.BDnloc[2],sou.BDnloc[3]-1] - sou.mt[6]*sou.waveform[it]

end

#=== 3 ===# # add single force source #
function addsf!(wf::Union{elwf3d,acwf3d}, sou::SFsource, isn::Int64, it::Int64)

    wf.vx[sou.BDnloc[1],sou.BDnloc[2],sou.BDnloc[3]] = wf.vx[sou.BDnloc[1],sou.BDnloc[2],sou.BDnloc[3]] + sou.coeff[1]*sou.waveform[it]

    wf.vy[sou.BDnloc[1],sou.BDnloc[2],sou.BDnloc[3]] = wf.vy[sou.BDnloc[1],sou.BDnloc[2],sou.BDnloc[3]] + sou.coeff[2]*sou.waveform[it]

    wf.vz[sou.BDnloc[1],sou.BDnloc[2],sou.BDnloc[3]] = wf.vz[sou.BDnloc[1],sou.BDnloc[2],sou.BDnloc[3]] + sou.coeff[3]*sou.waveform[it]

end
