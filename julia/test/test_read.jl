import Downloads

function read_VSPT()
    mkpath("test_files")
    if !isfile("test_files/VSPT_ASCII.xyz")
        Downloads.download(raw"https://nasa-public-data.s3.amazonaws.com/plot3d_utilities/VSPT_ASCII.xyz","test_files/VSPT_ASCII.xyz")
    end
        
    return "completed"
    # blocks = Block[]
    # if isfile(filename)
    #     io = open(filename,"r")
    #     IMAX = Int64[]; JMAX = Int64[]; KMAX = Int64[]
    #     NBLOCKS = Int64(read(io,Int32))
    #     for i in 1:NBLOCKS
    #         append!(IMAX,Int64(read(io,Int32)))
    #         append!(JMAX,Int64(read(io,Int32)))
    #         append!(KMAX,Int64(read(io,Int32)))
    #     end
    #     for i in 1:NBLOCKS
    #         X = read_plot3D_chunk(io,IMAX[i],JMAX[i],KMAX[i])
    #         Y = read_plot3D_chunk(io,IMAX[i],JMAX[i],KMAX[i])
    #         Z = read_plot3D_chunk(io,IMAX[i],JMAX[i],KMAX[i])
    #         b = Block(IMAX=IMAX[i],JMAX=JMAX[i],KMAX=KMAX[i],X=X,Y=Y,Z=Z)
    #         push!(blocks,b)
    #     end
    #     close(io)
    # end
    # return blocks
end