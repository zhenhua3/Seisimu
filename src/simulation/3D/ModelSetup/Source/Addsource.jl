function Addsource!(WF::WF2D, sn::Int64, source::Array{Source,1}, it::Int64)

    # it : total time samplings
    # sn : source number
    for isn = 1:sn
        WF.VecBDTxx[source[isn].BDposn[1,1]] = WF.VecBDTxx[source[isn].BDposn[1,1]] + source[isn].BDposn[2,1]*source[isn].waveform[it]
        WF.VecBDTzz[source[isn].BDposn[1,1]] = WF.VecBDTzz[source[isn].BDposn[1,1]] + source[isn].BDposn[2,1]*source[isn].waveform[it]
        WF.VecBDVx[source[isn].BDposn[1,2]] = WF.VecBDVx[source[isn].BDposn[1,2]] + source[isn].BDposn[3,2]*source[isn].waveform[it]
        WF.VecBDVx[source[isn].BDposn[1,3]] = WF.VecBDVx[source[isn].BDposn[1,3]] + source[isn].BDposn[3,3]*source[isn].waveform[it]
        WF.VecBDVz[source[isn].BDposn[1,4]] = WF.VecBDVz[source[isn].BDposn[1,4]] + source[isn].BDposn[4,4]*source[isn].waveform[it]
        WF.VecBDVz[source[isn].BDposn[1,5]] = WF.VecBDVz[source[isn].BDposn[1,5]] + source[isn].BDposn[4,5]*source[isn].waveform[it]
    end
end

function Addsource!(WF::WF2D, source::Source, it::Int64)

    # it : total time samplings
        WF.VecBDTxx[source.BDposn[1,1]] = WF.VecBDTxx[source.BDposn[1,1]] + source.BDposn[2,1]*source.waveform[it]
        WF.VecBDTzz[source.BDposn[1,1]] = WF.VecBDTzz[source.BDposn[1,1]] + source.BDposn[2,1]*source.waveform[it]
        WF.VecBDVx[source.BDposn[1,2]] = WF.VecBDVx[source.BDposn[1,2]] + source.BDposn[3,2]*source.waveform[it]
        WF.VecBDVx[source.BDposn[1,3]] = WF.VecBDVx[source.BDposn[1,3]] + source.BDposn[3,3]*source.waveform[it]
        WF.VecBDVz[source.BDposn[1,4]] = WF.VecBDVz[source.BDposn[1,4]] + source.BDposn[4,4]*source.waveform[it]
        WF.VecBDVz[source.BDposn[1,5]] = WF.VecBDVz[source.BDposn[1,5]] + source.BDposn[4,5]*source.waveform[it]
end

function AddSLsource!(WF::WF2D, sn::Int64, source::Array{Source,1}, it::Int64)

    # it : total time samplings
    # sn : source number
    for isn = 1:sn
        WF.VecBDTxxx[source[isn].BDposn[1,1]] = WF.VecBDTxxx[source[isn].BDposn[1,1]] + source[isn].BDposn[2,1]*source[isn].waveform[it]/2
        WF.VecBDTzzx[source[isn].BDposn[1,1]] = WF.VecBDTzzx[source[isn].BDposn[1,1]] + source[isn].BDposn[2,1]*source[isn].waveform[it]/2
        WF.VecBDVxx[source[isn].BDposn[1,2]] = WF.VecBDVxx[source[isn].BDposn[1,2]] + source[isn].BDposn[3,2]*source[isn].waveform[it]/2
        WF.VecBDVxx[source[isn].BDposn[1,3]] = WF.VecBDVxx[source[isn].BDposn[1,3]] + source[isn].BDposn[3,3]*source[isn].waveform[it]/2
        WF.VecBDVzx[source[isn].BDposn[1,4]] = WF.VecBDVzx[source[isn].BDposn[1,4]] + source[isn].BDposn[4,4]*source[isn].waveform[it]/2
        WF.VecBDVzx[source[isn].BDposn[1,5]] = WF.VecBDVzx[source[isn].BDposn[1,5]] + source[isn].BDposn[4,5]*source[isn].waveform[it]/2
        WF.VecBDTxxz[source[isn].BDposn[1,1]] = WF.VecBDTxxz[source[isn].BDposn[1,1]] + source[isn].BDposn[2,1]*source[isn].waveform[it]/2
        WF.VecBDTzzz[source[isn].BDposn[1,1]] = WF.VecBDTzzz[source[isn].BDposn[1,1]] + source[isn].BDposn[2,1]*source[isn].waveform[it]/2
        WF.VecBDVxz[source[isn].BDposn[1,2]] = WF.VecBDVxz[source[isn].BDposn[1,2]] + source[isn].BDposn[3,2]*source[isn].waveform[it]/2
        WF.VecBDVxz[source[isn].BDposn[1,3]] = WF.VecBDVxz[source[isn].BDposn[1,3]] + source[isn].BDposn[3,3]*source[isn].waveform[it]/2
        WF.VecBDVzz[source[isn].BDposn[1,4]] = WF.VecBDVzz[source[isn].BDposn[1,4]] + source[isn].BDposn[4,4]*source[isn].waveform[it]/2
        WF.VecBDVzz[source[isn].BDposn[1,5]] = WF.VecBDVzz[source[isn].BDposn[1,5]] + source[isn].BDposn[4,5]*source[isn].waveform[it]/2
    end
end

function AddSLsource!(WF::WF2D, source::Source, it::Int64)

    # it : total time samplings
    WF.VecBDTxxx[source.BDposn[1,1]] = WF.VecBDTxxx[source.BDposn[1,1]] + source.BDposn[2,1]*source.waveform[it]/2
    WF.VecBDTzzx[source.BDposn[1,1]] = WF.VecBDTzzx[source.BDposn[1,1]] + source.BDposn[2,1]*source.waveform[it]/2
    WF.VecBDVxx[source.BDposn[1,2]] = WF.VecBDVxx[source.BDposn[1,2]] + source.BDposn[3,2]*source.waveform[it]/2
    WF.VecBDVxx[source.BDposn[1,3]] = WF.VecBDVxx[source.BDposn[1,3]] + source.BDposn[3,3]*source.waveform[it]/2
    WF.VecBDVzx[source.BDposn[1,4]] = WF.VecBDVzx[source.BDposn[1,4]] + source.BDposn[4,4]*source.waveform[it]/2
    WF.VecBDVzx[source.BDposn[1,5]] = WF.VecBDVzx[source.BDposn[1,5]] + source.BDposn[4,5]*source.waveform[it]/2
    WF.VecBDTxxz[source.BDposn[1,1]] = WF.VecBDTxxz[source.BDposn[1,1]] + source.BDposn[2,1]*source.waveform[it]/2
    WF.VecBDTzzz[source.BDposn[1,1]] = WF.VecBDTzzz[source.BDposn[1,1]] + source.BDposn[2,1]*source.waveform[it]/2
    WF.VecBDVxz[source.BDposn[1,2]] = WF.VecBDVxz[source.BDposn[1,2]] + source.BDposn[3,2]*source.waveform[it]/2
    WF.VecBDVxz[source.BDposn[1,3]] = WF.VecBDVxz[source.BDposn[1,3]] + source.BDposn[3,3]*source.waveform[it]/2
    WF.VecBDVzz[source.BDposn[1,4]] = WF.VecBDVzz[source.BDposn[1,4]] + source.BDposn[4,4]*source.waveform[it]/2
    WF.VecBDVzz[source.BDposn[1,5]] = WF.VecBDVzz[source.BDposn[1,5]] + source.BDposn[4,5]*source.waveform[it]/2
end
