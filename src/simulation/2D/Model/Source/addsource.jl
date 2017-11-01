function addsou!(wf::nspwf2d, sou::Array{source}, it::Int64)
    sn = length(sou)
    # it : total time samplings
    # sn : sou number
        for isn = 1:sn
            wf.txx[sou[isn].BDnloc[1,1]] = wf.txx[sou[isn].BDnloc[1,1]] + sou[isn].BDnloc[2,1]*sou[isn].waveform[it]
            wf.tzz[sou[isn].BDnloc[1,1]] = wf.tzz[sou[isn].BDnloc[1,1]] + sou[isn].BDnloc[2,1]*sou[isn].waveform[it]
            wf.vx[sou[isn].BDnloc[1,2]] = wf.vx[sou[isn].BDnloc[1,2]] + sou[isn].BDnloc[3,2]*sou[isn].waveform[it]
            wf.vx[sou[isn].BDnloc[1,3]] = wf.vx[sou[isn].BDnloc[1,3]] + sou[isn].BDnloc[3,3]*sou[isn].waveform[it]
            wf.vz[sou[isn].BDnloc[1,4]] = wf.vz[sou[isn].BDnloc[1,4]] + sou[isn].BDnloc[4,4]*sou[isn].waveform[it]
            wf.vz[sou[isn].BDnloc[1,5]] = wf.vz[sou[isn].BDnloc[1,5]] + sou[isn].BDnloc[4,5]*sou[isn].waveform[it]
        end
end

function addsou!(wf::spwf2d, sou::Array{source,1}, it::Int64)
    sn = length(sou)
    # it : total time samplings
    # sn : sou number
        for isn = 1:sn
            wf.txxx[sou[isn].BDnloc[1,1]] = wf.txxx[sou[isn].BDnloc[1,1]] + sou[isn].BDnloc[2,1]*sou[isn].waveform[it]/2
            wf.tzzx[sou[isn].BDnloc[1,1]] = wf.tzzx[sou[isn].BDnloc[1,1]] + sou[isn].BDnloc[2,1]*sou[isn].waveform[it]/2
            wf.vxx[sou[isn].BDnloc[1,2]] = wf.vxx[sou[isn].BDnloc[1,2]] + sou[isn].BDnloc[3,2]*sou[isn].waveform[it]/2
            wf.vxx[sou[isn].BDnloc[1,3]] = wf.vxx[sou[isn].BDnloc[1,3]] + sou[isn].BDnloc[3,3]*sou[isn].waveform[it]/2
            wf.vzx[sou[isn].BDnloc[1,4]] = wf.vzx[sou[isn].BDnloc[1,4]] + sou[isn].BDnloc[4,4]*sou[isn].waveform[it]/2
            wf.vzx[sou[isn].BDnloc[1,5]] = wf.vzx[sou[isn].BDnloc[1,5]] + sou[isn].BDnloc[4,5]*sou[isn].waveform[it]/2
            wf.txxz[sou[isn].BDnloc[1,1]] = wf.txxz[sou[isn].BDnloc[1,1]] + sou[isn].BDnloc[2,1]*sou[isn].waveform[it]/2
            wf.tzzz[sou[isn].BDnloc[1,1]] = wf.tzzz[sou[isn].BDnloc[1,1]] + sou[isn].BDnloc[2,1]*sou[isn].waveform[it]/2
            wf.vxz[sou[isn].BDnloc[1,2]] = wf.vxz[sou[isn].BDnloc[1,2]] + sou[isn].BDnloc[3,2]*sou[isn].waveform[it]/2
            wf.vxz[sou[isn].BDnloc[1,3]] = wf.vxz[sou[isn].BDnloc[1,3]] + sou[isn].BDnloc[3,3]*sou[isn].waveform[it]/2
            wf.vzz[sou[isn].BDnloc[1,4]] = wf.vzz[sou[isn].BDnloc[1,4]] + sou[isn].BDnloc[4,4]*sou[isn].waveform[it]/2
            wf.vzz[sou[isn].BDnloc[1,5]] = wf.vzz[sou[isn].BDnloc[1,5]] + sou[isn].BDnloc[4,5]*sou[isn].waveform[it]/2
        end
end
