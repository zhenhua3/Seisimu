
# Default data precision is double

function TotalHToDDataSizeInBytes(nwf::nelwf3d)
    RowNumber = 0;
    vx_offset = nwf.BDnvx[1]*nwf.BDnvx[2]*nwf.BDnvx[3];
    RowNumber += vx_offset;
    vy_offset = nwf.BDnvy[1]*nwf.BDnvy[2]*nwf.BDnvy[3];
    RowNumber += vy_offset;
    vz_offset = nwf.BDnvz[1]*nwf.BDnvz[2]*nwf.BDnvz[3];
    RowNumber += vz_offset;
    tpp_offset = nwf.BDntpp[1]*nwf.BDntpp[2]*nwf.BDntpp[3];
    RowNumber += 6*tpp_offset; # including rho, lambda, mu, and three tpp
    txy_offset = nwf.BDntxy[1]*nwf.BDntxy[2]*nwf.BDntxy[3];
    RowNumber += txy_offset;
    txz_offset = nwf.BDntxz[1]*nwf.BDntxz[2]*nwf.BDntxz[3];
    RowNumber += txz_offset;
    tyz_offset = nwf.BDntyz[1]*nwf.BDntyz[2]*nwf.BDntyz[3];
    RowNumber += tyz_offset;

    return RowNumber*8
end


function TotalHToDDataSizeInBytes(nwf::nacwf3d)
     RowNumber = 0;
     vx_offset = nwf.BDnvx[1]*nwf.BDnvx[2]*nwf.BDnvx[3];
     RowNumber += vx_offset;
     vy_offset = nwf.BDnvy[1]*nwf.BDnvy[2]*nwf.BDnvy[3];
     RowNumber += vy_offset;
     vz_offset = nwf.BDnvz[1]*nwf.BDnvz[2]*nwf.BDnvz[3];
     RowNumber += vz_offset;
     tpp_offset = nwf.BDntpp[1]*nwf.BDntpp[2]*nwf.BDntpp[3];
     RowNumber += 3*tpp_offset;

     return RowNumber*8
end
