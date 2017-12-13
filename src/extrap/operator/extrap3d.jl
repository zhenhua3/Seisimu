function extrap3d!(model::elmod3d, nchunks::Int64)

    @sync for p in workers()
        @spawnat p begin
            vxbtxx!(model,nchunks)
            vxbtxz!(model,nchunks)
            vxbtxy!(model,nchunks)
            vzbtxz!(model,nchunks)
            vzbtzz!(model,nchunks)
            vzbtyz!(model,nchunks)
            vybtxy!(model,nchunks)
            vybtyy!(model,nchunks)
            vybtyz!(model,nchunks)
        end
    end
    @sync for p in workers()
        @spawnat p begin
            txxbvx!(model,nchunks)
            txxbvy!(model,nchunks)
            txxbvz!(model,nchunks)
            tyybvx!(model,nchunks)
            tyybvy!(model,nchunks)
            tyybvz!(model,nchunks)
            tzzbvx!(model,nchunks)
            tzzbvy!(model,nchunks)
            tzzbvz!(model,nchunks)
            txybvx!(model,nchunks)
            txybvy!(model,nchunks)
            txzbvx!(model,nchunks)
            txzbvz!(model,nchunks)
            tyzbvy!(model,nchunks)
            tyzbvz!(model,nchunks)
        end
    end
end


function extrap3d!(model::acmod3d, nchunks::Int64)

    @sync for p in workers()
        @spawnat p begin
            vxbtxx!(model,nchunks)
            vzbtzz!(model,nchunks)
            vybtyy!(model,nchunks)

        end
    end
    @sync for p in workers()
        @spawnat p begin
            txxbvx!(model,nchunks)
            txxbvy!(model,nchunks)
            txxbvz!(model,nchunks)
            tyybvx!(model,nchunks)
            tyybvy!(model,nchunks)
            tyybvz!(model,nchunks)
            tzzbvx!(model,nchunks)
            tzzbvy!(model,nchunks)
            tzzbvz!(model,nchunks)
        end
    end
end
